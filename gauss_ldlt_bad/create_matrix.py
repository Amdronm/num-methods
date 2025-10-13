import random
import subprocess
import time


def create_matrix(dim: int) -> list[list[float]]:
    mat: list[list[float]] = [[0] * dim for _ in range(dim)]
    for i in range(dim):
        for j in range(i, dim):
            mat[i][j] = random.uniform(-100_000.0, 100_000.0)
            mat[j][i] = mat[i][j]

    return mat


def write_matrix(mat: list[list[float]]):
    with open("input_rand.txt", "w") as f:
        f.write(str(len(mat)) + "\n")
        for row in mat:
            f.write(" ".join([str(elem) for elem in row]))
            f.write("\n")


def run_process() -> float:
    start = time.perf_counter()

    res = subprocess.run(args=["./build/gauss_ldlt", "input_rand.txt", "output.txt"])

    if res:
        end = time.perf_counter()
        return end - start
    else:
        return -1


def main():
    dim = 1000
    write_matrix(create_matrix(dim))
    run_process()


if __name__ == "__main__":
    main()
