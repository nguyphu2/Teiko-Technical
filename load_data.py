import argparse
import sqlite3
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import dash
from dash import dash_table, dcc, html, Input, Output
import plotly.express as px
import plotly.graph_objects as go


DEFAULT_DB = "cell_counts.db"

C = {
    "bg":     "#0f1923",
    "card":   "#162030",
    "border": "#1e3048",
    "accent": "#00c2ff",
    "yes":    "#00d4a8",
    "no":     "#ff6b6b",
    "text":   "#e8f4f8",
    "muted":  "#6b8fa8",
    "sig":    "#ffd166",
}

TABLE_STYLE = dict(
    style_table={"overflowX": "auto", "borderRadius": "8px", "border": f"1px solid {C['border']}"},
    style_header={"backgroundColor": C["border"], "color": C["accent"],
                  "fontWeight": "700", "fontSize": "12px", "letterSpacing": "0.08em",
                  "textTransform": "uppercase", "padding": "10px 14px", "border": "none"},
    style_cell={"backgroundColor": C["card"], "color": C["text"], "fontSize": "13px",
                "padding": "9px 14px", "border": f"1px solid {C['border']}",
                "fontFamily": "'IBM Plex Mono', monospace"},
    style_data_conditional=[
        {"if": {"row_index": "odd"}, "backgroundColor": "#192a3e"},
        {"if": {"filter_query": '{significant} = "True"'},
         "backgroundColor": "#2a2510", "color": C["sig"]},
    ],
)

PLOT_LAYOUT = dict(
    paper_bgcolor="rgba(0,0,0,0)",
    plot_bgcolor="rgba(0,0,0,0)",
    font=dict(family="'IBM Plex Mono', monospace", color=C["text"], size=12),
    xaxis=dict(gridcolor=C["border"], linecolor=C["border"], tickcolor=C["border"]),
    yaxis=dict(gridcolor=C["border"], linecolor=C["border"], tickcolor=C["border"]),
    legend=dict(bgcolor="rgba(0,0,0,0)", bordercolor=C["border"], borderwidth=1),
    margin=dict(t=50, b=50, l=60, r=20),
)


def get_connection(db_path):

    p = Path(db_path)

    if not p.exists():
        sys.exit(f"ERROR: database not found: {db_path}\nRun your main script first.")

    return sqlite3.connect(db_path)


def get_freq(conn):

    sql_query = """
        SELECT cc.sample_id AS sample,
               cc.cell_type AS population,
               cc.count     AS count
        FROM   cell_counts cc
        ORDER  BY cc.sample_id, cc.cell_type
    """

    df = pd.read_sql_query(sql_query, conn)

    totals = df.groupby("sample")["count"].sum().rename("total_count").reset_index()

    df = df.merge(totals, on="sample")

    df["percentage"] = (df["count"] / df["total_count"] * 100).round(2)

    return df[["sample", "total_count", "population", "count", "percentage"]]


def melanoma_pbmc(conn):

    sql_query = """
        SELECT s.sample_id                    AS sample,
               s.subject_id,
               sub.response,
               s.time_from_treatment_start,
               cc.cell_type                   AS population,
               cc.count                       AS count
        FROM   cell_counts cc
        JOIN   samples     s   ON cc.sample_id = s.sample_id
        JOIN   subjects    sub ON s.subject_id = sub.subject_id
        WHERE  sub.condition  = 'melanoma'
          AND  sub.treatment  = 'miraclib'
          AND  s.sample_type  = 'PBMC'
        ORDER  BY s.sample_id, cc.cell_type
    """

    df = pd.read_sql_query(sql_query, conn)

    totals = df.groupby("sample")["count"].sum().rename("total_count").reset_index()

    df = df.merge(totals, on="sample")

    df["percentage"] = (df["count"] / df["total_count"] * 100).round(2)

    return df


def FDR_correction(p_values):

    n = len(p_values)
    indexed = sorted(enumerate(p_values), key=lambda x: x[1])
    adjusted = np.zeros(n)
    prev = 1.0

    for rank, (orig_idx, p) in enumerate(reversed(indexed), start=1):
        adj = p * n / (n - rank + 1)
        adj = min(adj, prev)
        prev = adj
        adjusted[orig_idx] = adj

    return adjusted


def cal_statistics(df):

    results = []

    for pop, grp in df.groupby("population"):
        resp    = grp.loc[grp["response"] == "yes", "percentage"].values
        nonresp = grp.loc[grp["response"] == "no",  "percentage"].values

        if len(resp) < 2 or len(nonresp) < 2:
            u_stat, p_val = np.nan, np.nan
        else:
            u_stat, p_val = stats.mannwhitneyu(resp, nonresp, alternative="two-sided")

        results.append({
            "population":            pop,
            "n_responders":          len(resp),
            "n_non_responders":      len(nonresp),
            "mean_responders":       round(float(np.mean(resp)),    2) if len(resp)    else np.nan,
            "mean_non_responders":   round(float(np.mean(nonresp)), 2) if len(nonresp) else np.nan,
            "median_responders":     round(float(np.median(resp)),    2) if len(resp)    else np.nan,
            "median_non_responders": round(float(np.median(nonresp)), 2) if len(nonresp) else np.nan,
            "U_statistic":           round(u_stat, 3) if not np.isnan(u_stat) else np.nan,
            "p_value":               round(p_val,  4) if not np.isnan(p_val)  else np.nan,
        })

    rdf = pd.DataFrame(results)

    valid = rdf["p_value"].notna()
    p_adj = np.full(len(rdf), np.nan)
    p_adj[valid] = FDR_correction(rdf.loc[valid, "p_value"].tolist())
    rdf["p_adj"] = np.round(p_adj, 4)
    rdf["significant"] = rdf["p_adj"] < 0.05

    return rdf.sort_values("p_value").reset_index(drop=True)


def baseline_query(conn):

    sql_query = """
        SELECT s.sample_id,
               s.subject_id,
               sub.project,
               sub.condition,
               sub.treatment,
               sub.sex,
               sub.response,
               s.sample_type,
               s.time_from_treatment_start
        FROM   samples  s
        JOIN   subjects sub ON s.subject_id = sub.subject_id
        WHERE  sub.condition               = 'melanoma'
          AND  sub.treatment               = 'miraclib'
          AND  s.sample_type               = 'PBMC'
          AND  s.time_from_treatment_start = 0
        ORDER  BY sub.project, s.sample_id
    """

    return pd.read_sql_query(sql_query, conn)


def samples_per_project(df):
    return df.groupby("project").size().rename("sample_count").reset_index().sort_values("project")


def subjects_by_response(df):

    unique_subjects = df.drop_duplicates(subset="subject_id")

    counts = (
        unique_subjects.groupby("response")
        .size()
        .rename("subject_count")
        .reset_index()
    )

    counts["response"] = counts["response"].map({"yes": "Responder", "no": "Non-responder"})

    return counts.sort_values("response")


def subjects_by_sex(df):

    unique_subjects = df.drop_duplicates(subset="subject_id")

    return (
        unique_subjects.groupby("sex")
        .size()
        .rename("subject_count")
        .reset_index()
        .sort_values("sex")
    )


def load_all(db_path):

    conn = get_connection(db_path)

    freq_df  = get_freq(conn)
    mel_df   = melanoma_pbmc(conn)
    stats_df = cal_statistics(mel_df)
    baseline = baseline_query(conn)

    conn.close()

    return {
        "freq":        freq_df,
        "mel":         mel_df,
        "stats":       stats_df,
        "baseline":    baseline,
        "by_project":  samples_per_project(baseline),
        "by_response": subjects_by_response(baseline),
        "by_sex":      subjects_by_sex(baseline),
    }


def make_boxplot(mel_df):

    fig = px.box(
        mel_df,
        x="population", y="percentage", color="response",
        color_discrete_map={"yes": C["yes"], "no": C["no"]},
        labels={"population": "Cell Population", "percentage": "Relative Frequency (%)",
                "response": "Response"},
        category_orders={"population": sorted(mel_df["population"].unique())},
        points="all",
    )

    fig.update_traces(marker_size=6, line_width=1.5)

    fig.update_layout(
        **PLOT_LAYOUT,
        title=dict(text="Cell Population Frequencies — Responders vs Non-Responders",
                   font_size=14, x=0.01),
        height=440,
        boxgap=0.25,
    )

    return fig


def make_freq_bar(freq_df, selected_sample):

    sub = freq_df[freq_df["sample"] == selected_sample]

    fig = px.bar(
        sub, x="population", y="percentage", color="population",
        text=sub["percentage"].apply(lambda v: f"{v:.1f}%"),
        labels={"population": "Cell Type", "percentage": "Frequency (%)"},
        color_discrete_sequence=px.colors.qualitative.Bold,
    )

    fig.update_traces(textposition="outside", marker_line_width=0)

    fig.update_layout(
        **PLOT_LAYOUT,
        showlegend=False,
        title=dict(text=f"Frequency breakdown — {selected_sample}", font_size=13, x=0.01),
        height=360,
        yaxis_range=[0, sub["percentage"].max() * 1.2],
    )

    return fig


def make_pie(df, names_col, values_col, title):

    pie_layout = {
        **PLOT_LAYOUT,
        "margin": dict(t=45, b=20, l=20, r=20),
        "legend": dict(orientation="h", y=-0.1, bgcolor="rgba(0,0,0,0)",
                       bordercolor=C["border"], borderwidth=1),
    }

    fig = px.pie(
        df, names=names_col, values=values_col,
        color_discrete_sequence=[C["yes"], C["no"], C["accent"], C["sig"]],
    )

    fig.update_traces(textfont_size=13, marker=dict(line=dict(color=C["bg"], width=2)))

    fig.update_layout(**pie_layout, title=dict(text=title, font_size=13, x=0.01), height=300)

    return fig


def card(*children, style=None):

    base = {
        "backgroundColor": C["card"],
        "borderRadius": "10px",
        "border": f"1px solid {C['border']}",
        "padding": "20px",
        "marginBottom": "20px",
    }

    if style:
        base.update(style)

    return html.Div(children, style=base)


def section_title(text):

    return html.H3(text, style={
        "color": C["accent"],
        "fontFamily": "'IBM Plex Mono', monospace",
        "fontSize": "13px",
        "letterSpacing": "0.12em",
        "textTransform": "uppercase",
        "marginBottom": "14px",
        "borderBottom": f"1px solid {C['border']}",
        "paddingBottom": "8px",
    })


def stat_badge(label, value):

    return html.Div([
        html.Div(str(value), style={
            "fontSize": "28px",
            "fontWeight": "700",
            "color": C["accent"],
            "fontFamily": "'IBM Plex Mono', monospace",
        }),
        html.Div(label, style={
            "fontSize": "11px",
            "color": C["muted"],
            "marginTop": "2px",
            "textTransform": "uppercase",
            "letterSpacing": "0.08em",
        }),
    ], style={
        "backgroundColor": C["bg"],
        "borderRadius": "8px",
        "padding": "16px 20px",
        "border": f"1px solid {C['border']}",
        "textAlign": "center",
        "flex": "1",
    })


def build_app(data):

    freq_df  = data["freq"]
    mel_df   = data["mel"]
    stats_df = data["stats"]
    baseline = data["baseline"]

    BASE_FONT = "'IBM Plex Sans', sans-serif"
    MONO_FONT = "'IBM Plex Mono', monospace"

    app = dash.Dash(
        __name__,
        title="Loblaw Bio — Trial Dashboard",
        external_stylesheets=[
            "https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500;700&family=IBM+Plex+Sans:wght@300;400;600&display=swap"
        ],
    )

    samples = sorted(freq_df["sample"].unique())

    tab2 = html.Div([
        card(
            section_title("Overview — all samples"),
            dash_table.DataTable(
                id="freq-table",
                columns=[{"name": c, "id": c} for c in freq_df.columns],
                data=freq_df.to_dict("records"),
                page_size=15,
                sort_action="native",
                filter_action="native",
                **TABLE_STYLE,
            ),
        ),
        card(
            section_title("Per-sample breakdown"),
            html.Div([
                html.Label("Select sample:", style={
                    "color": C["muted"], "fontSize": "12px",
                    "fontFamily": MONO_FONT, "marginBottom": "6px",
                }),
                dcc.Dropdown(
                    id="sample-picker",
                    options=[{"label": s, "value": s} for s in samples],
                    value=samples[0],
                    clearable=False,
                    style={
                        "backgroundColor": C["bg"], "color": C["text"],
                        "border": f"1px solid {C['border']}", "borderRadius": "6px",
                        "fontFamily": MONO_FONT, "fontSize": "13px",
                    },
                ),
            ], style={"maxWidth": "340px", "marginBottom": "16px"}),
            dcc.Graph(id="freq-bar"),
        ),
    ])

    stats_display = stats_df.copy()
    stats_display["significant"] = stats_display["significant"].astype(str)

    sig_pops = stats_df[stats_df["significant"]]["population"].tolist()

    sig_note = (
        f"Significant populations (p_adj < 0.05): {', '.join(sig_pops)}"
        if sig_pops else
        "No populations reached significance at p_adj < 0.05 — may reflect small sample size."
    )

    tab3 = html.Div([
        card(
            section_title("Responders vs Non-Responders — boxplots"),
            html.P("Melanoma patients · miraclib · PBMC samples only", style={
                "color": C["muted"], "fontSize": "12px",
                "fontFamily": MONO_FONT, "marginBottom": "14px",
            }),
            dcc.Graph(figure=make_boxplot(mel_df), id="boxplot"),
        ),
        card(
            section_title("Statistical Results — Mann-Whitney U + BH-FDR"),
            html.Div(sig_note, style={
                "backgroundColor": C["bg"],
                "border": f"1px solid {C['sig']}",
                "borderRadius": "6px",
                "padding": "10px 14px",
                "color": C["sig"],
                "fontSize": "12px",
                "fontFamily": MONO_FONT,
                "marginBottom": "16px",
            }),
            dash_table.DataTable(
                columns=[{"name": c, "id": c} for c in stats_display.columns],
                data=stats_display.to_dict("records"),
                sort_action="native",
                **TABLE_STYLE,
            ),
        ),
    ])

    n_samples  = len(baseline)
    n_subjects = baseline["subject_id"].nunique()
    n_projects = baseline["project"].nunique()

    tab4 = html.Div([
        card(
            section_title("Baseline cohort — melanoma · miraclib · PBMC · timepoint 0"),
            html.Div([
                stat_badge("Baseline Samples", n_samples),
                stat_badge("Unique Subjects",  n_subjects),
                stat_badge("Projects",         n_projects),
            ], style={"display": "flex", "gap": "16px", "marginBottom": "0"}),
        ),
        html.Div([
            html.Div(card(
                section_title("Samples per project"),
                dcc.Graph(
                    figure=make_pie(data["by_project"], "project", "sample_count", ""),
                    config={"displayModeBar": False},
                ),
            ), style={"flex": "1"}),
            html.Div(card(
                section_title("Subjects by response"),
                dcc.Graph(
                    figure=make_pie(data["by_response"], "response", "subject_count", ""),
                    config={"displayModeBar": False},
                ),
            ), style={"flex": "1"}),
            html.Div(card(
                section_title("Subjects by sex"),
                dcc.Graph(
                    figure=make_pie(data["by_sex"], "sex", "subject_count", ""),
                    config={"displayModeBar": False},
                ),
            ), style={"flex": "1"}),
        ], style={"display": "flex", "gap": "16px"}),
        card(
            section_title("Baseline sample details"),
            dash_table.DataTable(
                columns=[{"name": c, "id": c} for c in baseline.columns],
                data=baseline.to_dict("records"),
                sort_action="native",
                filter_action="native",
                page_size=10,
                **TABLE_STYLE,
            ),
        ),
    ])

    app.layout = html.Div([
        html.Div([
            html.Div([
                html.Span("LOBLAW BIO", style={
                    "color": C["accent"], "fontFamily": MONO_FONT,
                    "fontSize": "11px", "letterSpacing": "0.2em", "fontWeight": "700",
                }),
                html.H1("Clinical Trial Dashboard", style={
                    "color": C["text"], "fontFamily": BASE_FONT,
                    "fontWeight": "300", "fontSize": "26px", "margin": "4px 0 2px",
                }),
                html.P("miraclib · immune cell population analysis", style={
                    "color": C["muted"], "fontFamily": MONO_FONT,
                    "fontSize": "12px", "margin": "0",
                }),
            ]),
        ], style={
            "padding": "28px 32px 20px",
            "borderBottom": f"1px solid {C['border']}",
            "marginBottom": "24px",
        }),
        html.Div([
            dcc.Tabs(
                id="tabs",
                value="tab2",
                children=[
                    dcc.Tab(label="Part 2 — Frequencies",        value="tab2"),
                    dcc.Tab(label="Part 3 — Responder Analysis", value="tab3"),
                    dcc.Tab(label="Part 4 — Baseline Subset",    value="tab4"),
                ],
                colors={"border": C["border"], "primary": C["accent"], "background": C["card"]},
                style={"fontFamily": MONO_FONT, "fontSize": "12px", "letterSpacing": "0.05em"},
            ),
            html.Div(id="tab-content", style={"padding": "24px 0 0"}),
        ], style={"padding": "0 32px 40px"}),
    ], style={
        "backgroundColor": C["bg"],
        "minHeight": "100vh",
        "color": C["text"],
        "fontFamily": BASE_FONT,
    })

    @app.callback(Output("tab-content", "children"), Input("tabs", "value"))
    def render_tab(tab):
        if tab == "tab2": return tab2
        if tab == "tab3": return tab3
        if tab == "tab4": return tab4

    @app.callback(Output("freq-bar", "figure"), Input("sample-picker", "value"))
    def update_bar(sample):
        return make_freq_bar(freq_df, sample)

    return app


def main():

    parser = argparse.ArgumentParser(description="Loblaw Bio — Clinical Trial Dashboard")
    parser.add_argument("--db",    default=DEFAULT_DB, help="SQLite database file")
    parser.add_argument("--port",  default=8050, type=int)
    parser.add_argument("--debug", action="store_true")
    args = parser.parse_args()

    print(f"\nLoading data from {args.db} ...")

    data = load_all(args.db)

    print(f"  Part 2: {len(data['freq'])} frequency rows")
    print(f"  Part 3: {len(data['mel'])} melanoma-miraclib-PBMC rows")
    print(f"  Part 4: {len(data['baseline'])} baseline samples")
    print(f"\nDashboard running at http://127.0.0.1:{args.port}/\n")

    app = build_app(data)
    app.run(debug=args.debug, port=args.port)


if __name__ == "__main__":
    main()