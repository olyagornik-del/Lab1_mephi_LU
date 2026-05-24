from experiments import (
    experiment_single_system,
    experiment_multiple_rhs,
    experiment_hilbert,
)


def main():
    print("Сравнение методов решения СЛАУ")
    print("Все эксперименты воспроизводимы: используется фиксированный seed = 42")

    # если компьютер не тянет большие размеры, можно уменьшить
    single_sizes = [100, 200, 500]
    multiple_rhs_size = 500
    multiple_rhs_counts = [1, 10, 100]
    hilbert_sizes = [5, 10, 15]

    experiment_single_system(single_sizes, seed=42)
    experiment_multiple_rhs(multiple_rhs_size, multiple_rhs_counts, seed=42)
    experiment_hilbert(hilbert_sizes)


if __name__ == "__main__":
    main()
