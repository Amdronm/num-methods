import matplotlib.pyplot as plt


def plot_convergence():
    ns = [10, 100, 1000]
    colors = ["navy", "crimson", "forestgreen"]
    markers = ["o", "s", "^"]

    plt.figure(figsize=(7, 5))

    for n, color, marker in zip(ns, colors, markers):
        filename = f"conv{n}.txt"

        with open(filename, "r") as f:
            res = [float(line.strip()) for line in f if line.strip()]

        k = range(1, len(res) + 1)
        plt.plot(
            k,
            res,
            marker=marker,
            color=color,
            linestyle="-",
            linewidth=1.5,
            markersize=6,
            label=f"n = {n}",
        )

    plt.yscale("log")
    plt.xlabel("k")
    plt.ylabel("$||F(x^k)||$")
    plt.title("Диаграмма сходимости метода Ньютона")
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend()
    plt.tight_layout()
    plt.savefig("task2_convergence.png", dpi=300)
    plt.close()


def plot_solutions():
    ns = [10, 100, 1000]
    colors = ["navy", "crimson", "forestgreen"]

    for n, color in zip(ns, colors):
        filename = f"sol{n}.txt"

        with open(filename, "r") as f:
            lines = [line.split() for line in f if line.strip()]
            i_vals = [int(line[0]) for line in lines]
            x_vals = [float(line[1]) for line in lines]

        plt.figure(figsize=(6, 4.5))

        if n == 10:
            plt.plot(
                i_vals,
                x_vals,
                marker="o",
                color=color,
                linestyle="-",
                markersize=6,
                linewidth=1.5,
                label=f"n = {n}",
            )
        elif n == 100:
            plt.plot(
                i_vals,
                x_vals,
                marker=".",
                color=color,
                linestyle="-",
                markersize=4,
                linewidth=1.5,
                label=f"n = {n}",
            )
        else:
            plt.plot(i_vals, x_vals, color=color, linestyle="-", linewidth=2, label=f"n = {n}")

        plt.xlabel("i")
        plt.ylabel("$x_i$")
        plt.title(f"Решение системы (n = {n})")
        plt.grid(True, ls="--", alpha=0.5)
        plt.legend()

        out_filename = f"task2_solution_{n}.png"
        plt.tight_layout()
        plt.savefig(out_filename, dpi=300)


if __name__ == "__main__":
    plot_convergence()
    plot_solutions()
