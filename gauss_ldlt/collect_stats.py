from typing import Tuple
from create_matrix import create_matrix, write_matrix, run_process

import matplotlib.pyplot as plt


def run_benches() -> list[Tuple[int, float]]:
    sizes: list[int] = [100 + i * 100 for i in range(20)]
    result: list[Tuple[int, float]] = []
    for size in sizes:
        write_matrix(create_matrix(size))
        time = run_process()
        result.append((size, time))

    return result


def build_plot(bench_results: list[Tuple[int, float]]):
    valid_sizes = []
    valid_times = []

    for size, time in bench_results:
        if time > 0:
            valid_sizes.append(size)
            valid_times.append(time)

    plt.figure(figsize=(12, 8))

    plt.subplot(2, 1, 1)
    plt.plot(valid_sizes, valid_times, "bo-", linewidth=2, markersize=6, label="Execution Time")
    plt.xlabel("Matrix Size (n x n)")
    plt.ylabel("Time (seconds)")
    plt.title("Matrix Inversion Performance: Time vs Size")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.savefig("pictures/benches.png", dpi=300, bbox_inches="tight")


def main():
    build_plot(run_benches())


if __name__ == "__main__":
    main()
