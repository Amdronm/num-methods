import os

import matplotlib.pyplot as plt
import numpy as np


def read_residuals_file(filename):
    """
    Чтение файла со списком float (норма невязки)
    Каждая строка содержит одно число
    """
    with open(filename, "r") as f:
        residuals = [float(line.strip()) for line in f if line.strip()]
    return residuals


def read_coordinates_file(filename):
    """
    Чтение файла с координатами
    В каждой строке: номер строки i, значение x_i
    """
    indices = []
    values = []
    with open(filename, "r") as f:
        for line_num, line in enumerate(f, 1):
            if line.strip():
                parts = line.strip().split()
                if len(parts) >= 2:
                    # Предполагаем формат: i x_i
                    indices.append(int(parts[0]))
                    values.append(float(parts[1]))
                elif len(parts) == 1:
                    # Если только значение без индекса, используем номер строки
                    indices.append(line_num)
                    values.append(float(parts[0]))
    return indices, values


def plot_combined_residuals(file_list, labels=None):
    """
    Построение совмещенного графика норм невязок из нескольких файлов
    """
    plt.figure(figsize=(12, 8))

    if labels is None:
        labels = [f"Файл {i + 1}" for i in range(len(file_list))]

    colors = ["blue", "red", "green", "orange", "purple", "brown"]
    markers = ["o", "s", "^", "D", "v", "<"]

    for i, (filename, label) in enumerate(zip(file_list, labels)):
        try:
            residuals = read_residuals_file(filename)
            iterations = range(1, len(residuals) + 1)

            color = colors[i % len(colors)]
            marker = markers[i % len(markers)]

            plt.plot(
                iterations,
                residuals,
                color=color,
                marker=marker,
                markersize=3,
                linewidth=1.5,
                label=f"{label} (n={len(residuals)})",
                alpha=0.8,
            )

            print(f"Загружен {filename}: {len(residuals)} значений")

        except Exception as e:
            print(f"Ошибка при чтении {filename}: {e}")

    plt.xlabel("Итерация", fontsize=12)
    plt.ylabel("Норма невязки", fontsize=12)
    plt.title("Сравнение норм невязок по итерациям", fontsize=14)
    plt.legend(loc="best", fontsize=10)
    plt.grid(True, alpha=0.3, linestyle="--")

    # Логарифмическая шкала для лучшей визуализации
    plt.yscale("log")
    plt.ylabel("Норма невязки (логарифмическая шкала)", fontsize=12)

    plt.tight_layout()
    plt.savefig("some2.png")


def plot_individual_coordinates(file_list, labels=None):
    """
    Построение отдельных графиков для каждого файла с координатами
    """
    if labels is None:
        labels = [f"Файл {i + 1}" for i in range(len(file_list))]

    n_files = len(file_list)
    # Создаем сетку подграфиков
    fig, axes = plt.subplots(1, n_files, figsize=(5 * n_files, 5))

    if n_files == 1:
        axes = [axes]

    colors = ["blue", "red", "green", "orange", "purple", "brown"]

    for i, (filename, label) in enumerate(zip(file_list, labels)):
        try:
            indices, values = read_coordinates_file(filename)

            ax = axes[i]
            ax.plot(
                indices,
                values,
                color=colors[i % len(colors)],
                marker="o",
                markersize=4,
                linewidth=1.5,
                alpha=0.7,
            )

            ax.set_xlabel("Индекс i", fontsize=10)
            ax.set_ylabel("Значение x_i", fontsize=10)
            ax.set_title(f"{label}\n(n={len(values)})", fontsize=11)
            ax.grid(True, alpha=0.3, linestyle="--")

            print(f"Загружен {filename}: {len(values)} точек")

        except Exception as e:
            print(f"Ошибка при чтении {filename}: {e}")
            axes[i].text(
                0.5,
                0.5,
                f"Ошибка: {e}",
                ha="center",
                va="center",
                transform=axes[i].transAxes,
            )

    plt.tight_layout()
    plt.savefig("some.png")


def plot_combined_coordinates(file_list, labels=None):
    """
    Построение совмещенного графика координат из нескольких файлов
    """
    plt.figure(figsize=(12, 8))

    if labels is None:
        labels = [f"Файл {i + 1}" for i in range(len(file_list))]

    colors = ["blue", "red", "green", "orange", "purple", "brown"]
    markers = ["o", "s", "^", "D", "v", "<"]

    for i, (filename, label) in enumerate(zip(file_list, labels)):
        try:
            indices, values = read_coordinates_file(filename)

            color = colors[i % len(colors)]
            marker = markers[i % len(markers)]

            plt.plot(
                indices,
                values,
                color=color,
                marker=marker,
                markersize=3,
                linewidth=1.5,
                label=f"{label} (n={len(values)})",
                alpha=0.7,
                linestyle="-",
            )

            print(f"Загружен {filename}: {len(values)} точек")

        except Exception as e:
            print(f"Ошибка при чтении {filename}: {e}")

    plt.xlabel("Индекс i", fontsize=12)
    plt.ylabel("Значение x_i", fontsize=12)
    plt.title("Сравнение координат из разных файлов", fontsize=14)
    plt.legend(loc="best", fontsize=10)
    plt.grid(True, alpha=0.3, linestyle="--")
    plt.tight_layout()
    plt.savefig("some.png")


def save_data_info(file_list, file_type="residuals"):
    """
    Вывод информации о загруженных файлах
    """
    print(f"\n{'=' * 50}")
    print(f"Информация о {file_type} файлах:")
    print(f"{'=' * 50}")

    for i, filename in enumerate(file_list):
        try:
            if file_type == "residuals":
                data = read_residuals_file(filename)
                stats = {
                    "min": np.min(data),
                    "max": np.max(data),
                    "mean": np.mean(data),
                    "std": np.std(data),
                    "first": data[0],
                    "last": data[-1],
                }
                print(f"\nФайл {i + 1}: {filename}")
                print(f"  Размер: {len(data)}")
                print(f"  Диапазон: [{stats['min']:.6e}, {stats['max']:.6e}]")
                print(f"  Среднее: {stats['mean']:.6e} ± {stats['std']:.6e}")
                print(f"  Первое/последнее: {stats['first']:.6e} / {stats['last']:.6e}")
            else:
                indices, values = read_coordinates_file(filename)
                stats = {
                    "min": np.min(values),
                    "max": np.max(values),
                    "mean": np.mean(values),
                    "std": np.std(values),
                }
                print(f"\nФайл {i + 1}: {filename}")
                print(f"  Размер: {len(values)}")
                print(f"  Диапазон x: [{stats['min']:.6e}, {stats['max']:.6e}]")
                print(f"  Среднее x: {stats['mean']:.6e} ± {stats['std']:.6e}")
                print(f"  Индексы: от {indices[0]} до {indices[-1]}")
        except Exception as e:
            print(f"  Ошибка: {e}")


# ============= ОСНОВНАЯ ПРОГРАММА =============

if __name__ == "__main__":
    # ===== ЧАСТЬ 1: Файлы с нормами невязок =====
    # Укажите имена ваших 3 файлов с нормами невязок
    residual_files = [
        "conv10.txt",  # Замените на реальные имена файлов
        "conv100.txt",
        "conv1000.txt",
    ]

    # Подписи для графиков (опционально)
    residual_labels = ["Метод 1", "Метод 2", "Метод 3"]

    # ===== ЧАСТЬ 2: Файлы с координатами =====
    # Укажите имена ваших 3 файлов с координатами
    coordinate_files = [
        "sol10.txt",  # Замените на реальные имена файлов
    ]

    # Подписи для графиков координат (опционально)
    coordinate_labels = ["Набор данных 1", "Набор данных 2", "Набор данных 3"]

    # Проверяем существование файлов
    print("Проверка файлов...")

    # Для невязок
    existing_residuals = [f for f in residual_files if os.path.exists(f)]
    if len(existing_residuals) < len(residual_files):
        print(
            f"Предупреждение: найдено {len(existing_residuals)} из {len(residual_files)} файлов с невязками"
        )
        residual_files = existing_residuals
        residual_labels = residual_labels[: len(residual_files)]

    if residual_files:
        print(f"\nНайдены файлы с невязками: {residual_files}")
        save_data_info(residual_files, "residuals")

        # Строим совмещенный график невязок
        plot_combined_residuals(residual_files, residual_labels)
    else:
        print("Не найдено ни одного файла с невязками!")

    # Для координат
    existing_coords = [f for f in coordinate_files if os.path.exists(f)]
    if len(existing_coords) < len(coordinate_files):
        print(
            f"\nПредупреждение: найдено {len(existing_coords)} из {len(coordinate_files)} файлов с координатами"
        )
        coordinate_files = existing_coords
        coordinate_labels = coordinate_labels[: len(coordinate_files)]

    if coordinate_files:
        print(f"\nНайдены файлы с координатами: {coordinate_files}")
        save_data_info(coordinate_files, "coordinates")

        # Строим отдельные графики для каждого файла координат
        plot_individual_coordinates(coordinate_files, coordinate_labels)

        # Также можно построить совмещенный график координат
        # plot_combined_coordinates(coordinate_files, coordinate_labels)
    else:
        print("Не найдено ни одного файла с координатами!")
