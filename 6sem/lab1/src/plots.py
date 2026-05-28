import os

import matplotlib.pyplot as plt
import pandas as pd


def plot_task1():
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    roots = [1, 2, 3]
    intervals = ["[0.1, 0.3]", "[0.3, 0.4]", "[0.4, 0.6]"]

    for i, root in enumerate(roots):
        file = f"task1_root{root}.csv"
        if os.path.exists(file):
            df = pd.read_csv(file)
            axes[i].plot(
                df["Iter"], df["Bisection_Res"], marker="o", label="Метод бисекции"
            )
            axes[i].plot(
                df["Iter"],
                df["Newton_Res"],
                marker="s",
                label="Метод Ньютона (const df)",
            )
            axes[i].set_yscale("log")
            axes[i].set_xlabel("Номер итерации $k$")
            axes[i].set_ylabel("Невязка $|f(x_k)|$")
            axes[i].set_title(f"Сходимость к корню {root} на {intervals[i]}")
            axes[i].grid(True, which="both", ls="--")
            axes[i].legend()
    plt.tight_layout()
    plt.savefig("task1_convergence.png", dpi=300)
    plt.show()


def plot_task2():
    ns = [10, 100, 1000]

    # 1. График сходимости для системы
    plt.figure(figsize=(8, 6))
    for n in ns:
        file = f"task2_conv_{n}.csv"
        if os.path.exists(file):
            df = pd.read_csv(file)
            plt.plot(df["Iter"], df["NormF"], marker="o", label=f"n = {n}")

    plt.yscale("log")
    plt.xlabel("Номер итерации $k$")
    plt.ylabel("Невязка $||f(x^k)||_\infty$")
    plt.title("Диаграмма сходимости метода Ньютона для СЛАУ")
    plt.grid(True, which="both", ls="--")
    plt.legend()
    plt.savefig("task2_convergence.png", dpi=300)
    plt.show()

    # 2. Графики решений (точечные)
    plt.figure(figsize=(8, 6))
    for n in ns:
        file = f"task2_sol_{n}.csv"
        if os.path.exists(file):
            df = pd.read_csv(file)
            plt.plot(
                df["i"] * (10.0 / n),
                df["xi"],
                label=f"n = {n}",
                markersize=3,
                linestyle="-" if n == 1000 else "None",
                marker="o" if n < 1000 else "None",
            )

    plt.xlabel("Координата $x$ (масштабировано до [0,10])")
    plt.ylabel("Значение $x_i$")
    plt.title("Решения СЛАУ при разных n")
    plt.grid(True)
    plt.legend()
    plt.savefig("task2_solutions.png", dpi=300)
    plt.show()


if __name__ == "__main__":
    plot_task1()
    plot_task2()
