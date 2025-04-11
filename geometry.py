import numpy as np
import plotly.graph_objects as go
from ipywidgets import (
    FloatText, Button, Dropdown, VBox, HBox, Output, IntText, Layout, Tab
)

def update_plot(option, wire, qpcgap, qpcheight, vdiag_val, nx, ny, n=None):
    delx, dely = 1.0, 1.0
    xmin, ymin = 1.0, 1.0
    xmax, ymax = delx * nx, dely * ny
    def vdiag_qpc(ix, iy):
        x, y = xmin + delx * ix, ymin + dely * iy
        if (ymax - qpcheight) / 2 <= y <= (ymax + qpcheight) / 2:
            if (xmax - qpcgap) / 2 <= x <= (xmax + qpcgap) / 2:
                return 0.0
        elif wire <= x <= xmax - wire:
            return 0.0
        return vdiag_val

    def vdiag_dot(ix, iy):
        x, y = xmin + delx * ix, ymin + dely * iy
        centers = [(i + 1) * ymax / (n + 2) for i in range(n + 1)] if n else []
        for center in centers:
            if ((center - qpcheight / 2 <= y <= center + qpcheight / 2) and
                (x <= (xmax - qpcgap) / 2 or x >= (xmax + qpcgap) / 2)):
                return vdiag_val
        if x <= wire or x >= xmax - wire:
            return vdiag_val
        return 0.0

    def vdiag_one_sided_dot(ix, iy):
        x, y = xmin + delx * ix, ymin + dely * iy
        centers = [(i + 1) * ymax / (n + 2) for i in range(n + 1)] if n else []
        for center in centers:
            if ((center - qpcheight / 2 <= y <= center + qpcheight / 2) and
                (wire + qpcgap <= x <= xmax - wire)):
                return vdiag_val
        if x <= wire or x >= xmax - wire:
            return vdiag_val
        return 0.0

    vdiag_func = {
        "Quantum Point Contact": vdiag_qpc,
        "Quantum Dot": vdiag_dot,
        "One-Sided Quantum Dot": vdiag_one_sided_dot
    }.get(option, vdiag_qpc)

    potential_vals = np.zeros((nx, ny))
    for ix in range(nx):
        for iy in range(ny):
            potential_vals[ix, iy] = vdiag_func(ix, iy)

    x_vals = np.linspace(xmin, xmax, nx)
    y_vals = np.linspace(ymin, ymax, ny)
    x_grid, y_grid = np.meshgrid(x_vals, y_vals, indexing='ij')

    with open("fort.45", "w") as f:
        for ix in range(nx):
            for iy in range(ny):
                f.write(f"{x_grid[ix, iy]:.6f} {y_grid[ix, iy]:.6f} {potential_vals[ix, iy]:.6f}\n")

    x_range = xmax - xmin
    y_range = ymax - ymin
    z_range = potential_vals.max() - potential_vals.min()

    aspect_ratio = dict(
        x=1,
        y=y_range / x_range,
        z=0.15*z_range
    )

    fig = go.Figure(data=[
        go.Surface(
            z=potential_vals,
            x=x_grid,
            y=y_grid,
            colorscale='Viridis',
            colorbar=dict(title="Potential")
        )
    ])

    fig.update_layout(
        title="3D Potential Distribution",
        scene=dict(
            xaxis_title="X-axis",
            yaxis_title="Y-axis",
            zaxis_title="Potential",
            aspectratio=aspect_ratio,
        ),
        width=700,
        height=600,
    )

    with output:
        output.clear_output(wait=True)
        fig.show()

option_input = Dropdown(
    options=["Quantum Point Contact", "Quantum Dot", "One-Sided Quantum Dot"],
    description="Type:",
    value="Quantum Point Contact"
)
wire_input = FloatText(value=10.0, description='Wire:', layout=Layout(width='150px'))
qpcgap_input = FloatText(value=10.0, description='QPC Gap:', layout=Layout(width='150px'))
qpcheight_input = FloatText(value=15.0, description='QPC Length:', layout=Layout(width='150px'))
vdiag_val_input = FloatText(value=3.0, description='Vdiag Value:', layout=Layout(width='150px'))
n_input = IntText(value=1, description='Num Dots:', layout=Layout(width='150px'))
n_input.layout.display = 'none'

nx_input = IntText(value=100, description="X Grid Size:", layout=Layout(width='150px'))
ny_input = IntText(value=200, description="Y Grid Size:", layout=Layout(width='150px'))

def on_option_change(change):
    n_input.layout.display = '' if change.new in ["Quantum Dot", "One-Sided Quantum Dot"] else 'none'

option_input.observe(on_option_change, names='value')

submit_button = Button(description="Generate Plot", button_style='primary')
output = Output()

def click(b):
    update_plot(
        option_input.value,
        wire_input.value,
        qpcgap_input.value,
        qpcheight_input.value,
        vdiag_val_input.value,
        nx_input.value,
        ny_input.value,
        n_input.value if option_input.value in ["Quantum Dot", "One-Sided Quantum Dot"] else None
    )

submit_button.on_click(click)

tabs = Tab(children=[
    VBox([HBox([option_input, wire_input, qpcgap_input, qpcheight_input, vdiag_val_input, n_input], layout=Layout(align_items='center'))]),
    VBox([HBox([nx_input, ny_input])])
])
tabs.set_title(0, "Geometry")
tabs.set_title(1, "Mesh Size")

layout = VBox([tabs, submit_button, output])
layout
