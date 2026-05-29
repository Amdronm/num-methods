import matplotlib.pyplot as plt


def read_results(filename):
    residuals = []

    with open(filename, "r") as f:
        lines = f.read().splitlines()

    for i in range(0, len(lines), 3):
        if i + 1 >= len(lines):
            break
        count = int(lines[i + 1].strip())
        res = list(map(float, lines[i + 2].split())) if count > 0 else []
        residuals.append(res)
    return residuals


def main():
    res_b = read_results("outputB.txt")
    res_n = read_results("outputN.txt")
    intervals = ["[0.2182736, 0.3]", "[0.3, 0.45]", "[0.45, 0.54789798]"]

    for i in range(3):
        plt.figure(figsize=(6, 4.5))

        if i < len(res_b) and res_b[i]:
            plt.plot(
                range(1, len(res_b[i]) + 1),
                res_b[i],
                "o-",
                color="navy",
                label="Бисекция",
            )

        if i < len(res_n) and res_n[i]:
            plt.plot(
                range(1, len(res_n[i]) + 1),
                res_n[i],
                "s-",
                color="crimson",
                label="Ньютон",
            )

        plt.yscale("log")
        plt.xlabel("Номер итерации k")
        plt.ylabel("Невязка |f(x_k)|")
        plt.title(f"Интервал {intervals[i]}")
        plt.grid(True, which="both", ls="--", alpha=0.5)
        plt.legend()

        filename = f"plot{i + 1}.png"
        plt.tight_layout()
        plt.savefig(filename, dpi=150)


if __name__ == "__main__":
    main()
