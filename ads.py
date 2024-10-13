import numpy as np

# Входная последовательность в бинарной форме (GF(2))
sequence = np.array([1, 1, 1, 0, 1, 0, 0])

# Алгоритм Берлекэмпа-Мэсси для нахождения минимального многочлена
def berlekamp_massey(s):
    n = len(s)
    c = np.zeros(n, dtype=int)  # текущий многочлен
    b = np.zeros(n, dtype=int)  # предыдущий многочлен
    c[0] = b[0] = 1  # начальные условия
    l, m, d = 0, -1, 0  # степень, задержка, дискриминант

    for i in range(n):
        d = s[i]  # дискриминант
        for j in range(1, l + 1):
            d ^= c[j] * s[i - j]  # подсчёт дискриминанта
        if d == 1:  # требуется обновление многочлена
            temp = c.copy()
            p = np.roll(b, i - m)  # сдвиг предыдущего многочлена
            c ^= p  # обновление текущего многочлена
            if l <= i // 2:
                l = i + 1 - l
                m = i
                b = temp
    return c[:l+1], l

# Выполняем алгоритм
min_poly, degree = berlekamp_massey(sequence)

print(min_poly, degree)  # минимальный многочлен и его степень
