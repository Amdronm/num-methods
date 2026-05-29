import matplotlib.pyplot as plt


def read_residuals(filename):
    with open(filename, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    # Первая строка - количество интервалов
    # num_intervals = int(lines[0])

    residuals_by_interval = []
    values_by_interval = []

    idx = 0
    for interval in range(3):
        # Значение интервала
        interval_value = float(lines[idx])
        idx += 1

        # Количество невязок
        num_residuals = int(lines[idx])
        idx += 1

        # Строка с невязками
        residuals = list(map(float, lines[idx].split()))
        idx += 1

        # Проверка количества невязок
        if len(residuals) != num_residuals:
            print(
                f"Предупреждение: ожидалось {num_residuals} невязок, получено {len(residuals)}"
            )

        residuals_by_interval.append(residuals)
        values_by_interval.append(interval_value)

    return residuals_by_interval, values_by_interval


def plot_residuals(residuals_by_interval, values_by_interval):
    """
    Построение графика невязок по итерациям
    """
    plt.figure(figsize=(12, 8))

    colors = ["blue", "red", "green"]
    markers = ["o", "s", "^"]

    for i, (residuals, interval_value) in enumerate(
        zip(residuals_by_interval, values_by_interval)
    ):
        iterations = range(1, len(residuals) + 1)

        plt.plot(
            iterations,
            residuals,
            marker=markers[i],
            color=colors[i],
            markersize=4,
            linewidth=1.5,
            label=f"Интервал {i + 1}: {interval_value:.6f}",
        )

        # Добавление точек для наглядности
        plt.scatter(iterations, residuals, color=colors[i], s=20, alpha=0.6)

    plt.xlabel("Итерация", fontsize=12)
    plt.ylabel("Невязка", fontsize=12)
    plt.title("График невязок по итерациям для трех интервалов", fontsize=14)
    plt.legend(loc="best", fontsize=10)
    plt.grid(True, alpha=0.3, linestyle="--")

    # Логарифмическая шкала по y, если значения сильно отличаются
    # Проверяем, есть ли нулевые или отрицательные значения
    all_residuals = [r for residuals in residuals_by_interval for r in residuals]
    if min(all_residuals) > 0:
        plt.yscale("log")
        plt.ylabel("Невязка (логарифмическая шкала)", fontsize=12)

    plt.tight_layout()
    plt.show()


def plot_separate(residuals_by_interval, values_by_interval):
    """
    Построение отдельных графиков для каждого интервала
    """
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    colors = ["blue", "red", "green"]

    for i, (residuals, interval_value) in enumerate(
        zip(residuals_by_interval, values_by_interval)
    ):
        iterations = range(1, len(residuals) + 1)

        axes[i].plot(
            iterations,
            residuals,
            color=colors[i],
            marker="o",
            markersize=4,
            linewidth=1.5,
        )
        axes[i].set_xlabel("Итерация", fontsize=10)
        axes[i].set_ylabel("Невязка", fontsize=10)
        axes[i].set_title(f"Интервал {i + 1}: {interval_value:.6f}", fontsize=11)
        axes[i].grid(True, alpha=0.3, linestyle="--")

        # Добавляем логарифмическую шкалу при необходимости
        if min(residuals) > 0:
            axes[i].set_yscale("log")

    plt.tight_layout()
    plt.savefig("plottt.png")


# Основная программа
if __name__ == "__main__":
    # Имя файла с данными
    filename = "outputB.txt"  # Измените на имя вашего файла

    try:
        # Чтение данных из файла
        residuals_by_interval, values_by_interval = read_residuals(filename)

        # Вывод информации о данных
        print(f"Прочитано {len(residuals_by_interval)} интервалов")
        for i, (residuals, value) in enumerate(
            zip(residuals_by_interval, values_by_interval)
        ):
            print(
                f"Интервал {i + 1}: значение={value:.10f}, количество невязок={len(residuals)}"
            )
            print(f"  Невязки: первые 5 значений: {residuals[:5]}")

        # Построение общего графика
        plot_separate(residuals_by_interval, values_by_interval)

        # Построение отдельных графиков (опционально)
        # plot_separate(residuals_by_interval, values_by_interval)

    except FileNotFoundError:
        print(f"Ошибка: файл '{filename}' не найден")
    except Exception as e:
        print(f"Ошибка при обработке файла: {e}")
        print("\nПроверьте формат файла. Ожидаемый формат:")
        print("3")
        print("0.249999999994")
        print("30")
        print("0.000439247536779 0.00076885771204 ...")
