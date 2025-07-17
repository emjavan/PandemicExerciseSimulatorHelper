from dash import Dash, html, dcc, Input, Output, State, dash_table
import plotly.graph_objs as go
import numpy as np
from pymongo import MongoClient
import pandas as pd
import uuid

# --- Initialize MongoDB ---
client = MongoClient("mongodb://localhost:27017/")
db = client["TestSIRCompare"]
results_col = db["sir_results"]
results_col.delete_many({})  # clear previous runs each app start

# --- Define discrete SIR model ---
def discrete_sir(beta, gamma, I0, N, tolerance=0.9, max_days=365):
    S = N - I0
    I = I0
    R = 0

    S_list = [S]
    I_list = [I]
    R_list = [R]
    t_list = [0]

    day = 0
    while I > tolerance and day < max_days:
        new_infections = beta * S * I / N
        new_recoveries = gamma * I

        S -= new_infections
        I += new_infections - new_recoveries
        R += new_recoveries

        day += 1

        S_list.append(S)
        I_list.append(I)
        R_list.append(R)
        t_list.append(day)

    return {
        "time": np.array(t_list),
        "S": np.array(S_list),
        "I": np.array(I_list),
        "R": np.array(R_list),
        "stats": {
            "total_infected": R,
            "percent_infected": R / N * 100,
            "peak_infected": max(I_list),
            "days_to_peak": t_list[np.argmax(I_list)],
            "end_of_epidemic": day
        }
    }

# --- Initialize Dash app ---
app = Dash(__name__)

# --- Layout ---
app.layout = html.Div([
    html.H1("Baseline vs Intervention SIR Comparison"),

    html.H2("Baseline Parameters"),
    html.Label("Beta (infection rate) "),
    dcc.Input(id="beta_base", type="number", value=0.3, step=0.01),

    html.Label(" Gamma (recovery rate) "),
    dcc.Input(id="gamma_base", type="number", value=0.1, step=0.01),

    html.Label(" Initial Infected "),
    dcc.Input(id="I0_base", type="number", value=1, step=1),

    html.Label(" Population Size "),
    dcc.Input(id="N_base", type="number", value=1000, step=10),

    html.H2("Intervention Parameters"),
    html.Label("Beta (infection rate) "),
    dcc.Input(id="beta_int", type="number", value=0.2, step=0.01),

    html.Br(),
    html.Button("Run Comparison", id="run-button", n_clicks=0),

    html.Br(),
    html.Button("Export All Runs Summary", id="export-all-button", n_clicks=0),
    dcc.Download(id="download-all"),

    html.Br(),
    html.Button("Export All Time Series", id="export-timeseries-button", n_clicks=0),
    dcc.Download(id="download-timeseries"),

    dcc.Graph(id="sir-plot"),

    html.H2("Summary Statistics"),
    html.Div(id="summary-table")
])

# --- Callback to run baseline + intervention simulations ---
@app.callback(
    Output("sir-plot", "figure"),
    Output("summary-table", "children"),
    Input("run-button", "n_clicks"),
    State("beta_base", "value"),
    State("gamma_base", "value"),
    State("I0_base", "value"),
    State("N_base", "value"),
    State("beta_int", "value")
)
def run_comparison(n_clicks, beta_base, gamma_base, I0_base, N_base, beta_int):
    if n_clicks == 0:
        return go.Figure(), ""

    run_id = str(uuid.uuid4())

    # --- Run baseline ---
    baseline = discrete_sir(beta_base, gamma_base, I0_base, N_base)

    # --- Run intervention ---
    intervention = discrete_sir(beta_int, gamma_base, I0_base, N_base)

    # --- Save both runs to MongoDB ---
    results_col.insert_one({
        "run_id": run_id,
        "scenario": "baseline",
        "parameters": {"beta": beta_base, "gamma": gamma_base, "I0": I0_base, "N": N_base},
        "summary_stats": baseline["stats"],
        "time_series": {
            "time": baseline["time"].tolist(),
            "S": baseline["S"].tolist(),
            "I": baseline["I"].tolist(),
            "R": baseline["R"].tolist()
        }
    })

    results_col.insert_one({
        "run_id": run_id,
        "scenario": "intervention",
        "parameters": {"beta": beta_int, "gamma": gamma_base, "I0": I0_base, "N": N_base},
        "summary_stats": intervention["stats"],
        "time_series": {
            "time": intervention["time"].tolist(),
            "S": intervention["S"].tolist(),
            "I": intervention["I"].tolist(),
            "R": intervention["R"].tolist()
        }
    })

    # --- Plot baseline (solid) vs intervention (dashed) ---
    colors = {"Susceptible": "blue", "Infected": "red", "Recovered": "green"}
    fig = go.Figure()

    for compartment, data in zip(["Susceptible", "Infected", "Recovered"],
                                 [baseline["S"], baseline["I"], baseline["R"]]):
        fig.add_trace(go.Scatter(
            x=baseline["time"], y=data,
            mode="lines",
            name=f"Baseline {compartment}",
            line=dict(dash="solid", color=colors[compartment])
        ))

    for compartment, data in zip(["Susceptible", "Infected", "Recovered"],
                                 [intervention["S"], intervention["I"], intervention["R"]]):
        fig.add_trace(go.Scatter(
            x=intervention["time"], y=data,
            mode="lines",
            name=f"Intervention {compartment}",
            line=dict(dash="dash", color=colors[compartment])
        ))

    fig.update_layout(
        title="Baseline vs Intervention SIR",
        xaxis_title="Time",
        yaxis_title="Population",
        plot_bgcolor="white",
        paper_bgcolor="white",
        xaxis=dict(showgrid=True, gridcolor='lightgrey'),
        yaxis=dict(showgrid=True, gridcolor='lightgrey')
    )

    # --- Build summary table with difference ---
    stats_base = baseline["stats"]
    stats_int = intervention["stats"]

    data_table = [
        {
            "Metric": "Total Infected",
            "Baseline": round(stats_base["total_infected"],2),
            "Intervention": round(stats_int["total_infected"],2),
            "Difference": round(stats_int["total_infected"] - stats_base["total_infected"],2),
            "Run ID": run_id
        },
        {
            "Metric": "Percent Infected",
            "Baseline": round(stats_base["percent_infected"],2),
            "Intervention": round(stats_int["percent_infected"],2),
            "Difference": round(stats_int["percent_infected"] - stats_base["percent_infected"],2),
            "Run ID": run_id
        },
        {
            "Metric": "Peak Infected",
            "Baseline": round(stats_base["peak_infected"],2),
            "Intervention": round(stats_int["peak_infected"],2),
            "Difference": round(stats_int["peak_infected"] - stats_base["peak_infected"],2),
            "Run ID": run_id
        },
        {
            "Metric": "Days to Peak",
            "Baseline": round(stats_base["days_to_peak"],2),
            "Intervention": round(stats_int["days_to_peak"],2),
            "Difference": round(stats_int["days_to_peak"] - stats_base["days_to_peak"],2),
            "Run ID": run_id
        },
        {
            "Metric": "End of Epidemic (days)",
            "Baseline": round(stats_base["end_of_epidemic"],2),
            "Intervention": round(stats_int["end_of_epidemic"],2),
            "Difference": round(stats_int["end_of_epidemic"] - stats_base["end_of_epidemic"],2),
            "Run ID": run_id
        }
    ]

    table = dash_table.DataTable(
        data=data_table,
        columns=[{"name": i, "id": i} for i in ["Metric", "Baseline", "Intervention", "Difference", "Run ID"]],
        style_table={'width': '90%'}
    )

    return fig, table

# --- Export all runs summary ---
@app.callback(
    Output("download-all", "data"),
    Input("export-all-button", "n_clicks"),
    prevent_initial_call=True
)
def export_all(n_clicks):
    cursor = results_col.find()
    data = []
    for doc in cursor:
        stats = doc.get("summary_stats", {})
        row = {
            "run_id": doc.get("run_id"),
            "scenario": doc["scenario"],
            "beta": doc["parameters"]["beta"],
            "gamma": doc["parameters"]["gamma"],
            "I0": doc["parameters"]["I0"],
            "N": doc["parameters"]["N"],
            "total_infected": stats.get("total_infected"),
            "percent_infected": stats.get("percent_infected"),
            "peak_infected": stats.get("peak_infected"),
            "days_to_peak": stats.get("days_to_peak"),
            "end_of_epidemic": stats.get("end_of_epidemic")
        }
        data.append(row)

    df = pd.DataFrame(data)
    return dcc.send_data_frame(df.to_csv, "sir_all_runs_summary.csv")

# --- Export all time series ---
@app.callback(
    Output("download-timeseries", "data"),
    Input("export-timeseries-button", "n_clicks"),
    prevent_initial_call=True
)
def export_timeseries(n_clicks):
    cursor = results_col.find()
    data = []
    for doc in cursor:
        ts = doc.get("time_series", {})
        run_id = doc.get("run_id")
        scenario = doc.get("scenario")
        for t, S, I, R in zip(ts.get("time", []), ts.get("S", []), ts.get("I", []), ts.get("R", [])):
            data.append({
                "run_id": run_id,
                "scenario": scenario,
                "time": t,
                "S": S,
                "I": I,
                "R": R
            })

    df = pd.DataFrame(data)
    return dcc.send_data_frame(df.to_csv, "sir_all_timeseries.csv")

# --- Run app ---
if __name__ == "__main__":
    app.run(debug=True)





