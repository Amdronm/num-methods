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
