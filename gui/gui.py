#!/usr/bin/env python3
"""Small GUI to load an inc_ns input file, run it, and watch convergence live.

Usage: python3 gui.py
"""

import os
import subprocess
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext, ttk

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INC_EXEC = os.path.join(PROJECT_ROOT, "build", "inc_exec")
INPUT_DIR = os.path.join(PROJECT_ROOT, "input")
POLL_MS = 300


def convergence_csv_for(input_path):
    """Mirror cxx_main.C's output_prefix logic: <base>.i -> output/<base>_out_convergence.csv"""
    base = os.path.splitext(os.path.basename(input_path))[0]
    return os.path.join(PROJECT_ROOT, "output", base + "_out_convergence.csv")


class RunnerApp:
    def __init__(self, root):
        self.root = root
        root.title("inc_ns runner")
        root.geometry("1000x600")

        self.input_path = None
        self.proc = None
        self.log_path = None
        self.csv_path = None
        self.csv_pos = 0
        self.time_steps = []
        self.u_errors = []
        self.v_errors = []
        self.t_errors = []

        self._build_menu()
        self._build_layout()

    # ---- UI construction ----

    def _build_menu(self):
        menubar = tk.Menu(self.root)
        file_menu = tk.Menu(menubar, tearoff=0)
        file_menu.add_command(label="Open...", command=self.on_open)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.root.quit)
        menubar.add_cascade(label="File", menu=file_menu)
        self.root.config(menu=menubar)

    def _build_layout(self):
        left = ttk.Frame(self.root)
        left.pack(side=tk.LEFT, fill=tk.BOTH, expand=False, padx=6, pady=6)

        controls = ttk.Frame(left)
        controls.pack(side=tk.TOP, fill=tk.X)

        self.run_button = ttk.Button(controls, text="Run", command=self.on_run, state=tk.DISABLED)
        self.run_button.pack(side=tk.LEFT)

        self.status_var = tk.StringVar(value="No input file loaded")
        ttk.Label(controls, textvariable=self.status_var).pack(side=tk.LEFT, padx=8)

        self.text = scrolledtext.ScrolledText(left, width=50, height=35, state=tk.DISABLED)
        self.text.pack(side=tk.TOP, fill=tk.BOTH, expand=True, pady=(6, 0))

        right = ttk.Frame(self.root)
        right.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=6, pady=6)

        self.figure = Figure(figsize=(6, 5))
        self.ax = self.figure.add_subplot(111)
        self.ax.set_xlabel("time step")
        self.ax.set_ylabel("L2 error")
        self.ax.set_yscale("log")
        (self.u_line,) = self.ax.plot([], [], marker="o", label="u L2 error")
        (self.v_line,) = self.ax.plot([], [], marker="s", label="v L2 error")
        (self.t_line,) = self.ax.plot([], [], marker="^", label="T L2 error")
        self.ax.legend()

        self.canvas = FigureCanvasTkAgg(self.figure, master=right)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    # ---- File -> Open ----

    def on_open(self):
        path = filedialog.askopenfilename(
            title="Open input file",
            initialdir=INPUT_DIR if os.path.isdir(INPUT_DIR) else PROJECT_ROOT,
            filetypes=[("inc_ns input files", "*.i"), ("All files", "*")],
        )
        if not path:
            return

        with open(path, "r") as f:
            content = f.read()

        self.text.config(state=tk.NORMAL)
        self.text.delete("1.0", tk.END)
        self.text.insert(tk.END, content)
        self.text.config(state=tk.DISABLED)

        self.input_path = path
        self.run_button.config(state=tk.NORMAL)
        self.status_var.set("Loaded: " + os.path.basename(path))

    # ---- Run ----

    def on_run(self):
        if self.proc is not None:
            return  # a run is already in progress
        if not self.input_path:
            messagebox.showerror("No input file", "Open an input file first.")
            return

        self.csv_path = convergence_csv_for(self.input_path)
        if os.path.exists(self.csv_path):
            os.remove(self.csv_path)
        self.csv_pos = 0
        self.time_steps = []
        self.u_errors = []
        self.v_errors = []
        self.t_errors = []
        self._update_plot()

        base = os.path.splitext(os.path.basename(self.input_path))[0]
        self.log_path = os.path.join(PROJECT_ROOT, "output", base + "_out_gui.log")
        os.makedirs(os.path.dirname(self.log_path), exist_ok=True)
        log_file = open(self.log_path, "w")

        self.proc = subprocess.Popen(
            [INC_EXEC, "-i", self.input_path],
            cwd=PROJECT_ROOT,
            stdout=log_file,
            stderr=subprocess.STDOUT,
        )
        self._log_file_handle = log_file

        self.run_button.config(state=tk.DISABLED)
        self.status_var.set("Running...")
        self.root.after(POLL_MS, self._poll)

    def _poll(self):
        self._read_new_csv_rows()

        rc = self.proc.poll()
        if rc is None:
            self.root.after(POLL_MS, self._poll)
            return

        # Process finished: one last read to catch anything written just before exit.
        self._read_new_csv_rows()
        self._log_file_handle.close()
        self.proc = None
        self.run_button.config(state=tk.NORMAL)

        if rc == 0:
            self.status_var.set("Finished (rc=0)")
        else:
            self.status_var.set("Failed (rc=%d)" % rc)
            self._show_log_tail()

    def _read_new_csv_rows(self):
        if not self.csv_path or not os.path.exists(self.csv_path):
            return

        with open(self.csv_path, "r") as f:
            f.seek(self.csv_pos)
            new_text = f.read()
            self.csv_pos = f.tell()

        if not new_text:
            return

        changed = False
        for line in new_text.splitlines():
            line = line.strip()
            if not line or line.startswith("time_step"):
                continue
            parts = line.split(",")
            # The solver writes time_step,u,v,T -- the T column arrived with the energy
            # equation. Accept any row with at least the three older fields rather than
            # an exact count, so a future column can't silently blank the plot again.
            if len(parts) < 3:
                continue
            try:
                step = int(parts[0])
                u_err, v_err = float(parts[1]), float(parts[2])
                t_err = float(parts[3]) if len(parts) > 3 else 0.0
            except ValueError:
                continue
            self.time_steps.append(step)
            self.u_errors.append(u_err)
            self.v_errors.append(v_err)
            self.t_errors.append(t_err)
            changed = True

        if changed:
            self._update_plot()

    @staticmethod
    def _log_safe(values):
        """Blank out non-positive values for the log y-axis.

        A run with no temperature B.C. anywhere holds T at -T_init, so its "error" is an
        exact 0 on every step; u/v can also hit exact 0 once fully converged. Those can't
        be drawn on a log scale and would distort relim(), whereas NaN is simply skipped.
        """
        return [v if v > 0.0 else float("nan") for v in values]

    def _update_plot(self):
        self.u_line.set_data(self.time_steps, self._log_safe(self.u_errors))
        self.v_line.set_data(self.time_steps, self._log_safe(self.v_errors))
        self.t_line.set_data(self.time_steps, self._log_safe(self.t_errors))
        self.ax.relim()
        self.ax.autoscale_view()
        self.canvas.draw_idle()

    def _show_log_tail(self, n_lines=30):
        try:
            with open(self.log_path, "r") as f:
                lines = f.readlines()
        except OSError:
            return
        tail = "".join(lines[-n_lines:])
        messagebox.showerror("Simulation failed", "Last %d lines of output:\n\n%s" % (n_lines, tail))


def main():
    root = tk.Tk()
    RunnerApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
