import os
import pandas as pd
import time
import json
from pathlib import Path
from shiny import App, render, ui, reactive

from modules.utils import extract_zip_to_tmp, list_replicates_from_outputs, build_dfs_by_replicate, reorganize_extracted_fastqs, list_treatments_from_outputs
from modules.runner import pair_fastqs, run_crispresso_parallel
from modules import plots

# ===================== UI =====================
app_ui = ui.page_fluid(
    ui.h2("CRISPResso2 Analyzer"),

    ui.layout_columns(
        ui.card(
            ui.card_header("1) Upload FASTQs (zip) & set run parameters"),

            # ZIP input grande con texto de ayuda
            ui.div(
                ui.input_file(
                    "zip_fastqs",
                    "FASTQs Files:",
                    accept=[".zip"],
                    multiple=False,
                    placeholder="Use the browse button or drag and drop a file here",
                    width="60%"
                ),
            ),

            ui.layout_columns(
                ui.div(
                    ui.input_text("amplicon_seq", "Amplicon Sequence:", placeholder="Paste the amplicon sequence here..."),
                    ui.input_text("coding_seq", "Coding Sequence:", placeholder="Paste the coding sequence here..."),
                    ui.input_text("guide_seq", "Guide Sequence:", placeholder="Paste the guide sequence here..."),
                ),
                ui.div(
                    ui.input_numeric("min_aln_score", "Minimum alignment score:", value=60),
                    ui.input_numeric("plot_window", "Plot window size:", value=20),
                ),
                col_widths=[4, 4],
            ),

            ui.input_action_button(
                "run_crispresso",
                "Run CRISPResso on uploaded FASTQs",
                class_="btn-primary",
            ),
            ui.hr(),
            ui.output_text_verbatim("run_log", placeholder=True),
        ),
        ui.card(
            ui.card_header("2) Analyze CRISPResso Results"),
            ui.input_checkbox_group("sel_replicas", "Select replicates:", choices=[]),
            ui.input_text_area("mut_seq", "Enter the MUT sequence:", rows=3, placeholder="Paste the MUT sequence here...", width="60%"),
            ui.input_text_area("wt_seq", "Enter the WT* sequence:", rows=3, placeholder="Paste the WT* sequence here...", width="60%"),
            ui.input_select("ref_treatment", "Reference treatment for sensitivity:", choices=[]),
            ui.output_ui("order_editor"),
        ),
        col_widths=[6, 6]
    ),

    ui.layout_columns(
        ui.card(
            ui.card_header("Results table"),
            ui.output_ui("download_ui"),
            ui.output_ui("tabla_dataframe"),
        ),
        col_widths=[12]
    ),

    ui.layout_columns(
        ui.layout_columns(
            ui.card(
                ui.card_header("Frameshift% Plot"),
                ui.output_plot("plot_frameshift"),
            ),
            ui.card(
                ui.card_header("Indels% Plot"),
                ui.output_plot("plot_indels"),
            ),
            col_widths=[12, 12],
        ),
        ui.card(
            ui.card_header("Sensitivity Plot"),
            ui.output_plot("plot_sensitivity"),
        ),
        col_widths=[6, 6]
    ),
)

# =================== Server ===================
def server(input, output, session):
    # Holders
    tmp_zip_dir = reactive.Value(None)    # keep temp dir alive
    tmp_zip_path = reactive.Value(None)   # remember extracted tmp path
    outputs_root = reactive.Value(None)   # remember extracted root path
    replicas_holder = reactive.Value([])  # remember replica names
    dfs_by_rep = reactive.Value({})       # {rep: df}
    run_messages = reactive.Value("")     # log messages from CRISPResso run
    #outputs_tmpref = reactive.Value(None)  # <-- mantiene vivo el TemporaryDirectory

    def _append_log(msg:str):
        current = run_messages.get() or ""
        # Append newline only if something is already there, then add the new message
        run_messages.set(current + ("\n" if current else "") + msg) 

    # Check inputs validity (cleaned sequences, non-empty and valid)
    def _clean(seq):
        return (seq or "").replace("\n", "").replace(" ", "").upper().strip()

    def _valid(seq):
        if not seq:
            return False
        allowed = set("ACGTN")  
        return set(seq).issubset(allowed)
    
    # Store cleaned sequences for use in multiple places
    @reactive.calc
    def seqs():
        wt = (input.wt_seq() or "").strip().upper()
        mut = (input.mut_seq() or "").strip().upper()
        return mut, wt
    
    # ========== Step 1: ZIP upload only (no run yet) ==========
    @reactive.effect
    @reactive.event(input.zip_fastqs)
    def _on_zip_uploaded():
        fileinfo = input.zip_fastqs()
        if not fileinfo:
            tmp_zip_dir.set(None)
            tmp_zip_path.set(None)
            _append_log("No ZIP uploaded.")
            return
        
        datapath = fileinfo[0]["datapath"]
        tmp_dir = extract_zip_to_tmp(datapath)
        tmp_zip_dir.set(tmp_dir)
        reorganize_extracted_fastqs(Path(tmp_dir.name))
        tmp_zip_path.set(Path(tmp_dir.name))

        _append_log(f"Uploaded ZIP: {fileinfo[0]['name']} extracted to temporary directory.")


    # ========== Step 2: Run CRISPResso2 when button is clicked ==========
    outputs_tmpref = reactive.Value(None)
    @reactive.effect
    @reactive.event(input.run_crispresso)
    def _run_crispresso():
        _append_log("Starting CRISPResso2 run...")

        root_path = tmp_zip_path.get()
        if not root_path:
            _append_log("No ZIP uploaded yet. Please upload a ZIP file first.")
            return
        
        amplicon_seq = _clean(input.amplicon_seq())
        coding_seq  = _clean(input.coding_seq())
        guide_seq    = _clean(input.guide_seq())
        min_aln = input.min_aln_score() or 60
        plot_win = input.plot_window() or 20

        _append_log(f"Amplicon length: {len(amplicon_seq)}")
        _append_log(f"Coding seq length: {len(coding_seq)}")
        _append_log(f"Guide seq length: {len(guide_seq)}")
        _append_log(f"Minimum alignment score: {min_aln}")
        _append_log(f"Plot window size: {plot_win}")

        if not _valid(amplicon_seq) or not _valid(coding_seq) or not _valid(guide_seq):
            _append_log("Invalid sequences. Please check your input.")
            return
        
        # Pair FASTQs and run CRISPResso
        pairs = pair_fastqs(root_path)
        if not pairs:
            ui.notification_show("No valid FASTQ pairs found in the uploaded ZIP.", type="error")
            return

        out_root = Path("/outputs") / f"run-{int(time.time())}"
        out_root.mkdir(parents=True, exist_ok=True)
        outputs_root.set(out_root)
        _append_log(f"Running CRISPResso2 into: {out_root}")

        try:
            for line in run_crispresso_parallel(
                pairs=pairs,
                amplicon=amplicon_seq,
                guide=guide_seq,
                coding_seq=coding_seq,
                out_root=out_root,
                default_min_aln_score=int(min_aln),
                plot_window_size=int(plot_win),
            ):
                if line:
                    _append_log(line.rstrip())
        
        except Exception as e:
            _append_log(f"CRISPResso2 run failed: {e!r}")
   
        _append_log("CRISPResso2 finished. Discovering replicates...")
        reps = list_replicates_from_outputs(out_root)
        replicas_holder.set(reps)
        ui.update_checkbox_group("sel_replicas", choices=reps, selected=reps)
        treats = list_treatments_from_outputs(out_root)
        ui.update_select("ref_treatment", choices=treats, selected=(treats[0] if treats else None))
        _append_log(f"Replicates found: {', '.join(reps) if reps else '(none)'}")
        # ui.notification_show("CRISPResso2 finished.", type="message", duration=4)


    # ========== Step 3: Build DataFrames when inputs ready ==========
    @reactive.calc
    def ready_to_build():
        # returns TRUE if we have outputs root, replicas, and valid sequences
        return bool(outputs_root.get() and len(replicas_holder.get() or []) > 0 and all(_valid(s) for s in seqs()))

    @reactive.effect
    def _build_dfs():
        if not ready_to_build():
            return
       #  _append_log("Ready to build DataFrames for analysis...")
        root = outputs_root.get()       # directory with CRISPResso output file
        reps = replicas_holder.get()    # list of replicate folders to process
        wt  = (input.wt_seq() or "").strip().upper()
        mut = (input.mut_seq() or "").strip().upper()               # mut and wt* sequences
        amplicon = (input.amplicon_seq() or "").strip().upper() 
        ref_treat = input.ref_treatment()
        dfs = build_dfs_by_replicate(root, reps, wt, mut, amplicon, None, ref_treat)
        # Save the resulting dictionary of dfs into a reactive value
        dfs_by_rep.set(dfs)             

    @reactive.calc
    def df_selected():
        dfs = dfs_by_rep.get() or {}
        selected = input.sel_replicas() or list(dfs.keys())  # <- auto-selección
        frames = [dfs[r] for r in selected if r in dfs and not dfs[r].empty]
        if not frames:
            return pd.DataFrame()

        out = pd.concat(frames, ignore_index=True)

        if "Replicate" not in out.columns and "Sample" in out.columns:
            out.insert(0, "Replicate", out["Sample"])

        if "Replicate" in out.columns:
            out = out.loc[:, ["Replicate"] + [c for c in out.columns if c != "Replicate"]]

        return out

    def _current_treatment_order(base):
        if not base:
            return []
        order = []
        try:
            dnd = input.order_dnd()
        except Exception:
            dnd = None
        if dnd:
            try:
                arr = json.loads(dnd) if isinstance(dnd, str) else dnd
            except Exception:
                arr = dnd.split(",") if isinstance(dnd, str) else []
            for t in arr or []:
                if t in base and t not in order:
                    order.append(t)
        if not order:
            for i in range(len(base)):
                eid = f"order_pos_{i+1}"
                try:
                    val = getattr(input, eid)()
                except Exception:
                    val = None
                if val and val in base and val not in order:
                    order.append(val)
        for t in base:
            if t not in order:
                order.append(t)
        return order
    
    # ========== Outputs ==========
    @output
    @render.text
    def run_log():
        return run_messages.get() or "No log messages yet."

    @output
    @render.ui
    def tabla_dataframe():
        # thresholds for highlighting
        threshold = {
            "Total Reads": 1000,
            "MUT%": 1.00,
            "Indel%": 20,
        }

        d = df_selected()
        if d.empty:
            return ui.HTML("<div class='alert alert-info'>No data</div>")
        d_fmt = d.copy()

        # Order treatments and replicates
        order_base = list(pd.unique(d_fmt["Treatment"])) if "Treatment" in d_fmt.columns else []
        order = _current_treatment_order(order_base) if order_base else []
        if "Treatment" in d_fmt.columns and order:
            d_fmt["Treatment"] = pd.Categorical(d_fmt["Treatment"], categories=order, ordered=True)
        if "Replicate" in d_fmt.columns:
            reps_order = sorted(pd.unique(d_fmt["Replicate"]).tolist())
            d_fmt["Replicate"] = pd.Categorical(d_fmt["Replicate"], categories=reps_order, ordered=True)
        if "Replicate" in d_fmt.columns and "Treatment" in d_fmt.columns:
            d_fmt = d_fmt.sort_values(["Replicate", "Treatment"], kind="mergesort")
        
        # Columns to format: % in the name and numeric columns
        percent_cols = [c for c in d_fmt.columns if "%" in str(c)]
        center_cols = set(percent_cols)
        for c in d_fmt.columns:
            if pd.api.types.is_numeric_dtype(d_fmt[c]):
                center_cols.add(c)

        def format_cell(val, col):
            if pd.isna(val):
                return ""
            text = val

            # add % sign for percent columns
            if col in percent_cols and isinstance(val, (int, float)):
                text = f"{val}%"

            # decide if we highlight
            highlight_class = ""
            th = threshold.get(col, None)
            if th is not None:
                try:
                    v_float = float(val)
                    if v_float < th:
                        highlight_class = " cell-low"
                except (TypeError, ValueError):
                    pass

            # center if needed
            if col in center_cols:
                return f"<span class='cell-center{highlight_class}'>{text}</span>"
            else:
                return f"<span class='{highlight_class.strip()}'>{text}</span>"

        # Apply formatting column-wise
        for c in d_fmt.columns:
            d_fmt[c] = d_fmt[c].map(lambda x, col=c: format_cell(x, col))

        # Build CSS style and HTML table
        style = """
        <style>
            .table thead th { text-align:center; }
            .cell-center { display:block; text-align:center; }
            .cell-low { background-color: #f8d7da; }  /* light red */
        </style>
        """
        html = d_fmt.to_html(
            escape=False,
            index=False,
            classes=["table","table-striped","table-hover","table-sm"]
        )
        return ui.HTML(style + html)
    
    @output
    @render.ui
    def download_ui():
        d = df_selected()
        if d.empty:
            return None
        return ui.download_button("download", "Download CSV", class_="btn-primary btn-sm")
    
    @session.download(
        id="download",
        filename="crispresso2_analysis_results.csv")
    def download():
        d = df_selected()
        if d.empty:
            return None
        yield d.to_csv(index=False)

    @output
    @render.ui
    def order_editor():
        d = df_selected()
        if d.empty or "Treatment" not in d.columns:
            return ui.HTML("<div class='alert alert-warning'>No treatments to order yet.</div>")
        base = list(pd.unique(d["Treatment"]))
        try:
            dnd_raw = input.order_dnd()
        except Exception:
            dnd_raw = None
        user_order = []
        if dnd_raw:
            try:
                user_order = json.loads(dnd_raw) if isinstance(dnd_raw, str) else dnd_raw
            except Exception:
                user_order = []
        render_order = [t for t in user_order if t in base] + [t for t in base if t not in user_order]
        style = "<style>.order-title{font-size:1rem;font-weight:500;margin:4px 0 8px}.dnd-row{display:grid;grid-template-columns:26px 260px;gap:8px;align-items:center;margin-bottom:6px}.order-num{text-align:right;font-weight:600}.dnd-item{border:1px solid #ced4da;border-radius:4px;padding:6px 8px;background:#fff}.dnd-row.dragging{opacity:.6}</style>"
        items_html = "".join([f"<div class='dnd-row' draggable='true' data-value='{t}'><span class='order-num'>{i+1}.</span><div class='dnd-item'><span class='dnd-label'>{t}</span></div></div>" for i, t in enumerate(render_order)])
        script = "<script>(function(){var wrap=document.getElementById('dnd-treatments');if(!wrap)return;var last=wrap.getAttribute('data-last')||'';function setOrder(order){var json=JSON.stringify(order);if(json===last)return;last=json;wrap.setAttribute('data-last',json);if(window.Shiny&&window.Shiny.setInputValue){window.Shiny.setInputValue('order_dnd', json, {priority:'event'});}}function updateNumbers(){var rows=[].slice.call(wrap.querySelectorAll('.dnd-row'));rows.forEach(function(el,idx){var num=el.querySelector('.order-num');if(num){num.textContent=(idx+1)+'.';}});var order=rows.map(function(el){return el.getAttribute('data-value');});setOrder(order);}var rows=[].slice.call(wrap.querySelectorAll('.dnd-row'));rows.forEach(function(r){r.addEventListener('dragstart',function(e){r.classList.add('dragging');e.dataTransfer.effectAllowed='move';});r.addEventListener('dragend',function(){r.classList.remove('dragging');updateNumbers();});});wrap.addEventListener('dragover',function(e){e.preventDefault();var dragging=wrap.querySelector('.dnd-row.dragging');if(!dragging)return;var after=getAfter(wrap,e.clientY);if(after==null){wrap.appendChild(dragging);}else{wrap.insertBefore(dragging,after);}});function getAfter(container,y){var els=[].slice.call(container.querySelectorAll('.dnd-row:not(.dragging)'));var closest={offset:-Infinity,el:null};els.forEach(function(child){var box=child.getBoundingClientRect();var offset=y-box.top-box.height/2;if(offset<0 && offset>closest.offset){closest={offset:offset,el:child};}});return closest.el;}})();</script>"
        return ui.TagList(ui.HTML(style), ui.p("Order of treatment:", class_="order-title"), ui.HTML(f"<div id='dnd-treatments'>{items_html}</div>"), ui.HTML(script))

    @reactive.effect
    @reactive.event(input.view_report)
    def _on_view_report():
        url = input.view_report()
        if url:
            _append_log(f"Opening report: {url}")
    
    @output
    @render.plot
    def plot_frameshift():
        d = df_selected()
        if d.empty:
            return None
        order_base = list(pd.unique(d["Treatment"])) if "Treatment" in d.columns else []
        order = _current_treatment_order(order_base) if order_base else []
        return plots.make_frameshift_plot(d, treat_order=order)
    
    @output
    @render.plot
    def plot_indels():
        d = df_selected()
        if d.empty:
            return None
        order_base = list(pd.unique(d["Treatment"])) if "Treatment" in d.columns else []
        order = _current_treatment_order(order_base) if order_base else []
        return plots.make_indels_plot(d, treat_order=order)

    @output
    @render.plot
    def plot_sensitivity():
        d = df_selected()
        if d.empty:
            return None
        order_base = list(pd.unique(d["Treatment"])) if "Treatment" in d.columns else []
        order = _current_treatment_order(order_base) if order_base else []
        return plots.make_sensitivity_plot(d, treat_order=order)
        

app = App(app_ui, server)
