Свинтилов Андрей Павлович, 4 группа
# Лабораторная работа 2: "Интерполяция"
## Задание

Написать программу, которая для данного набора точек $(x_i, y_i)$ строит среднеквадратичное приближение $\varphi$ по произвольному количеству базисных функций, указанных в варианте. Итоговую линейную задачу наименьших квадратов следует решать методом QR-разложения.

Провести вычислительный эксперимент.
- Построить логарифмическую диаграмму изменения среднеквадратичного отклонения $\frac{\sqrt{\sum\limits_{i=0}^{n}(\varphi(x_i)-y_i)^2}}{n+1}$ при изменении числа $n$ от $0$ до $100$ с шагом $5$. 
- Построить графики полученных приближений для $n=1,5,20, 50$. На графиках должны быть отмечены точки, по которым строится приближение.
**Составить отчет**, который содержит
- Ответы на вопросы:
	1. как ведет себя величина среднеквадратичного отклонения с ростом $n$?
	2. Какое минимальное отклонение может быть достигнуто с помощью вашей программы? Почему?
- Графики и диаграммы
- Код программы

#### Вариант Е
Смещенный тригонометрический базис Фурье на отрезке $[a, b] : \{1, \cos(n\pi x/L), \sin (n\pi x/L)\},$ где $L=\frac{b-a}{2},\ n=1,2,3,\dots$  Метод QR-разложения.

## Решение
#### Немного теории
Ищем $\varphi$ в форме $\varphi(x)=\sum\limits_{i=0}^{n}\alpha_{i}\varphi_{i}(x)$ 
$$
\begin{cases}
\varphi_{0}(x)= 1 \\
\varphi_{2k-1}(x)=\cos\left(\frac{n\pi x}{L}\right) \\
\varphi_{2k}(x)=\sin\left(\frac{n\pi x}{L}\right)
\end{cases}
,\quad k=1,\dots, n
$$
Обозначим за $\Phi$ матрицу из столбцов $u_{j}=[\varphi_{j}(x_{0}), \varphi_{j}(x_{1}),\dots,\varphi_{j}(x_{n})]^{T}$, $y=[y_{0},y_{1},\dots,y_{N}]$ 
Применяя метод QR-разложения:

![[Pasted image 20260529161930.png]]

получаем искомый вектор $\alpha$.
#### Код
Выданные данные - файл `data.csv`.
Реализовали алгоритм на `Python`:
```python
import numpy as np
import pandas as pd


def Interpolate(x, y, n):
    N = len(x) - 1
    m = 2 * n + 1
    a, b = np.min(x), np.max(x)
    L = (b - a) / 2.0

    if m > N + 1:
        return None, None

    phi = np.zeros((N + 1, m))
    phi[:, 0] = 1.0

    for k in range(1, n + 1):
        phi[:, 2 * k - 1] = np.cos(k * np.pi * x / L)
        phi[:, 2 * k] = np.sin(k * np.pi * x / L)

    B = np.hstack((phi, y.reshape(-1, 1)))
    _, R = np.linalg.qr(B)

    if R.shape[0] > m:
        u = R[m:, -1]
        norm_u = np.linalg.norm(u)
    else:
        norm_u = 0.0

    mse = norm_u / (N + 1)
    R_hat = R[:m, :m]
    y_hat = R[:m, -1]
    alpha = np.linalg.solve(R_hat, y_hat)

    return mse, alpha


def main():
    df = pd.read_csv("data.csv")
    x = df["x"].values
    y = df["y"].values

    a, b = np.min(x), np.max(x)
    L = (b - a) / 2.0
    x_fine = np.linspace(a, b, 500)

    with open("error_diag.txt", "w") as f_err:
        n_vals = sorted(list([1] + list(range(0, 105, 5))))
        for n in n_vals:
            mse, alpha = Interpolate(x, y, n)

            if mse is None:
                break

            f_err.write(f"{n} {mse:.12e}\n")

            if n in [1, 5, 20, 50]:
                y_fine = np.full_like(x_fine, alpha[0])
                for k in range(1, n + 1):
                    y_fine += alpha[2 * k - 1] * np.cos(k * np.pi * x_fine / L)
                    y_fine += alpha[2 * k] * np.sin(k * np.pi * x_fine / L)

                with open(f"interp_n_{n}.txt", "w") as f_approx:
                    for x_val, y_val in zip(x_fine, y_fine):
                        f_approx.write(f"{x_val:.6f} {y_val:.6f}\n")


if __name__ == "__main__":
    main()


```
Строим логарифмическую диаграмму изменения среднеквадратического отклонения и  графики полученных приближений:
```python
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

```

#### Графики

Получили требуемые графики

![[error_log.png|662]]

![[interp_1.png|665]]

![[interp_5.png|662]]

![[interp_20.png|664]]

![[interp_50.png|665]]

#### Ответы на вопросы:

1. как ведет себя величина среднеквадратичного отклонения с ростом $n$?
	Видим, что отклонение строго монотонно уменьшается. До $n \approx 15$ оно стремительно падает. После увеличение количества базисных функций не дает такого результата - базиса из 31 функции хватает, чтобы достаточно точно описать кривизну исследуемой $\varphi$.
2. Какое минимальное отклонение может быть достигнуто с помощью вашей программы? Почему?
	Теоретически для любого набора данных можно получить отклонение равное нулю (не считая погрешность округления) - если взять количество базисных функций $m$ равное количеству точек $N+1$ - получим задачу интерполяции, имеющее точное решение. Используемый метод QR разложения устойчивый, поэтому увеличение $n$ не приведет к увеличению погрешности.
