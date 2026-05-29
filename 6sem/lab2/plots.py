import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot_error_diagram():
    data = np.loadtxt("error_diag.txt")
    n_vals = data[:, 0]
    mse_vals = data[:, 1]

    plt.figure(figsize=(7, 5))
    plt.plot(
        n_vals,
        mse_vals,
        marker="o",
        color="navy",
        linestyle="-",
        linewidth=1.5,
    )

    plt.yscale("log")
    plt.ylim(top=1.5)

    y_ticks = [0.04, 0.06, 0.1, 0.2, 0.4, 0.6, 1.0]
    plt.yticks(y_ticks, [str(val) for val in y_ticks])

    plt.xlabel("n", fontsize=11)
    plt.ylabel("Среднеквадратичное отклонение", fontsize=11)
    plt.title("Диаграмма изменения среднеквадратичного отклонения", fontsize=12)
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig("error_log.png", dpi=300)
    plt.close()


def plot_approximations():
    df = pd.read_csv("data.csv")
    x_pts = df["x"].values
    y_pts = df["y"].values

    ns = [1, 5, 20, 50]
    colors = ["crimson", "orange", "forestgreen", "purple"]

    for n, color in zip(ns, colors):
        filename = f"interp_n_{n}.txt"

        approx_data = np.loadtxt(filename)
        x_line = approx_data[:, 0]
        y_line = approx_data[:, 1]

        plt.figure(figsize=(7, 5))

        plt.scatter(x_pts, y_pts, color="navy", s=15, alpha=0.5, label="данные")

        plt.plot(
            x_line,
            y_line,
            color=color,
            linewidth=2,
            label=f"приближение (n = {n})",
        )

        plt.xlabel("x", fontsize=11)
        plt.ylabel("y", fontsize=11)
        plt.title(f"Приближение МНК (n = {n})", fontsize=12)
        plt.grid(True, ls="--", alpha=0.5)
        plt.legend()
        plt.tight_layout()

        out_filename = f"interp_{n}.png"
        plt.savefig(out_filename, dpi=300)
        plt.close()


if __name__ == "__main__":
    plot_error_diagram()
    plot_approximations()


# import os

# import matplotlib.pyplot as plt
# import numpy as np
# import pandas as pd


# def plot_error_diagram():
#     if not os.path.exists("error_diag.txt"):
#         print("Файл error_diag.txt не найден! Сначала запустите вычислительный скрипт.")
#         return

#     data = np.loadtxt("error_diag.txt")
#     n_vals = data[:, 0]
#     mse_vals = data[:, 1]

#     # --- ИСПРАВЛЕНИЕ БАГА ЛОГАРИФМА ---
#     # Отфильтровываем значения, близкие к нулю (точную интерполяцию),
#     # чтобы линия не улетала вертикально вниз в минус бесконечность
#     # mask = mse_vals > 1e-10
#     # n_vals = n_vals[mask]
#     # mse_vals = mse_vals[mask]
#     # ----------------------------------

#     plt.figure(figsize=(7, 5))
#     plt.plot(
#         n_vals,
#         mse_vals,
#         marker="o",
#         color="navy",
#         linestyle="-",
#         linewidth=1.5,
#         label="Среднеквадратичное отклонение",
#     )

#     # Делаем логарифмический масштаб по оси Y
#     plt.yscale("log")

#     # Приподнимаем верхнюю границу, чтобы влезла отметка 1.0 (10^0)
#     plt.ylim(top=1.5)

#     # Явно задаем красивые отметки по оси Y (включая 1.0)
#     y_ticks = [0.04, 0.06, 0.1, 0.2, 0.4, 0.6, 1.0]
#     plt.yticks(y_ticks, [str(val) for val in y_ticks])

#     # Настраиваем засечки по оси X (добавили 1)
#     x_ticks = [0, 1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95]
#     plt.xticks(x_ticks, fontsize=10)

#     plt.xlabel("n", fontsize=11)
#     plt.ylabel("Среднеквадратичное отклонение", fontsize=11)
#     plt.title("Диаграмма изменения среднеквадратичного отклонения", fontsize=12)
#     plt.grid(True, which="both", ls="--", alpha=0.5)
#     plt.legend()
#     plt.tight_layout()
#     plt.savefig("error_log.png", dpi=300)
#     plt.close()
#     print("Сохранена исправленная диаграмма ошибок: error_log.png")


# def plot_approximations():
#     df = pd.read_csv("data.csv")
#     x_pts = df["x"].values
#     y_pts = df["y"].values

#     ns = [1, 5, 20, 50]
#     colors = ["crimson", "orange", "forestgreen", "purple"]

#     for n, color in zip(ns, colors):
#         filename = f"interp_n_{n}.txt"
#         if not os.path.exists(filename):
#             print(f"Файл {filename} не найден!")
#             continue

#         approx_data = np.loadtxt(filename)
#         x_line = approx_data[:, 0]
#         y_line = approx_data[:, 1]

#         plt.figure(figsize=(7, 5))

#         plt.scatter(x_pts, y_pts, color="navy", s=15, alpha=0.5, label="данные")

#         plt.plot(
#             x_line,
#             y_line,
#             color=color,
#             linewidth=2,
#             label=f"приближение (n = {n})",
#         )

#         plt.xlabel("x", fontsize=11)
#         plt.ylabel("y", fontsize=11)
#         plt.title(f"Приближение МНК (n = {n})", fontsize=12)
#         plt.grid(True, ls="--", alpha=0.5)
#         plt.legend()
#         plt.tight_layout()

#         out_filename = f"interp_{n}.png"
#         plt.savefig(out_filename, dpi=300)
#         plt.close()
#         print(f"Сохранен график приближения: {out_filename}")


# if __name__ == "__main__":
#     plot_error_diagram()
#     plot_approximations()
