# By Yun Losson 2025 summer
# Input : pO2 measurements and x and y coordinates
# Output: reference point and pO2 measurements and its coordinates within defined polygon

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.widgets import PolygonSelector, Button
from matplotlib.colors import LinearSegmentedColormap
import tkinter as tk
from tkinter import filedialog
import matplotlib
matplotlib.use('TkAgg')

def main():
    fig, ax = plt.subplots(figsize=(12, 9))
    plt.subplots_adjust(bottom=0.2)

    data = None
    polygon_selector = None
    polygon_path = None
    polygon_vertices = []
    center_point = None
    center_coords = None
    center_pO2 = None
    polygon_started = False

    cmap = LinearSegmentedColormap.from_list('pO2_cmap', ['blue', 'cyan', 'green', 'yellow', 'red'])

    def upload_csv():
        nonlocal data, polygon_selector, polygon_vertices, center_point, center_coords, center_pO2, polygon_path, polygon_started

        root = tk.Tk()
        root.withdraw()
        file_path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])

        if file_path:
            try:
                data = pd.read_csv(file_path, encoding='ISO-8859-1')
                if len(data.columns) < 3:
                    print("Error: CSV needs at least 3 columns (x, y, pO2)")
                    return

                if polygon_selector:
                    polygon_selector.disconnect_events()
                polygon_vertices.clear()
                polygon_path = None
                polygon_started = False
                if center_point:
                    center_point.remove()
                    center_point = None
                center_coords = None
                center_pO2 = None

                ax.clear()
                x = data.iloc[:, 0]
                y = data.iloc[:, 1]
                pO2 = data.iloc[:, 2]

                grid_size = 10
                x_rounded = (x / grid_size).round() * grid_size
                y_rounded = (y / grid_size).round() * grid_size

                data_rounded = pd.DataFrame({'x': x_rounded, 'y': y_rounded, 'pO2': pO2})
                pivot = data_rounded.pivot_table(index='y', columns='x', values='pO2', aggfunc='mean')
                pivot = pivot.sort_index(ascending=False)

                X, Y = np.meshgrid(pivot.columns, pivot.index)
                mesh = ax.pcolormesh(X, Y, pivot.values, cmap=cmap, shading='auto')
                plt.colorbar(mesh, ax=ax, label='pO₂ Level')

                ax.scatter(x, y, c=pO2, cmap=cmap, s=0.1, alpha=0.0)

                ax.set_xlabel('x (µm)')
                ax.set_ylabel('y (µm)')
                ax.set_title('Click Center Point First, Then Draw Polygon')
                ax.invert_yaxis()

                fig.canvas.mpl_connect('button_press_event', on_click_set_center)

                plt.draw()

            except Exception as e:
                print(f"Error loading file: {e}")

    def on_click_set_center(event):
        nonlocal center_point, center_coords, center_pO2, data, polygon_selector, polygon_started

        if data is None or event.inaxes != ax:
            return

        if center_coords is not None:
            return  # Don't allow redefining center after it's selected

        if event.button == 1:
            x_click, y_click = event.xdata, event.ydata
            distances = np.sqrt((data.iloc[:, 0] - x_click)**2 + (data.iloc[:, 1] - y_click)**2)
            nearest_idx = distances.idxmin()

            center_coords = (data.iloc[nearest_idx, 0], data.iloc[nearest_idx, 1])
            center_pO2 = data.iloc[nearest_idx, 2]

            center_point = ax.scatter(
                *center_coords, c='magenta', s=150, marker='*',
                edgecolor='black', linewidth=0.5,
                label=f'Center (pO₂: {center_pO2:.1f})'
            )
            ax.legend()
            plt.draw()

            # Now allow polygon drawing
            polygon_selector = PolygonSelector(ax, onselect_polygon)
            polygon_started = True

    def onselect_polygon(vertices):
        nonlocal polygon_vertices, polygon_path

        polygon_vertices = vertices
        if len(vertices) >= 3:
            polygon_path = Path(vertices)
            plt.draw()

    def extract_data():
        nonlocal data, polygon_path, center_coords, center_pO2

        if data is None:
            print("No data loaded")
            return
        if polygon_path is None or len(polygon_vertices) < 3:
            print("No valid polygon selected")
            return
        if center_coords is None:
            print("No center point selected")
            return

        points = data.iloc[:, :2].values
        inside = polygon_path.contains_points(points)
        selected_data = data[inside].copy()

        if len(selected_data) == 0:
            print("No points found within the selected polygon")
            return

        center_header = pd.DataFrame({
            data.columns[0]: [center_coords[0]],
            data.columns[1]: [center_coords[1]],
            data.columns[2]: [center_pO2],
            'Note': ['Center Point']
        })

        output_data = pd.concat([center_header, selected_data], ignore_index=True)

        root = tk.Tk()
        root.withdraw()
        save_path = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv")],
            initialfile="pO2_selected_data.csv"
        )

        if save_path:
            try:
                output_data.to_csv(save_path, index=False)
                print(f"Data saved to {save_path}")
                print(f"Center point saved at top: {center_coords}, pO₂: {center_pO2}")
            except Exception as e:
                print(f"Error saving file: {e}")

    ax_upload = plt.axes([0.25, 0.05, 0.2, 0.075])
    ax_extract = plt.axes([0.55, 0.05, 0.2, 0.075])

    btn_upload = Button(ax_upload, 'Upload CSV')
    btn_upload.on_clicked(lambda event: upload_csv())

    btn_extract = Button(ax_extract, 'Extract Polygon Data')
    btn_extract.on_clicked(lambda event: extract_data())

    plt.show()

if __name__ == "__main__":
    main()
