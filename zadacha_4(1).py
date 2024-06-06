import numpy as np
import sympy as sym
import xlsxwriter


dannie = [[1, 0], [2, -1], [3, -3], [4, 2], [5, -2]]
eps = 0.000001
delta = eps / 4
workbook = xlsxwriter.Workbook('КР4_1.xlsx')
x0 = 10
y0 = 10


def f1(a, b):
    global dannie
    otvet = 0
    for i in range(len(dannie)):
        x = dannie[i][0]
        y = dannie[i][1]
        otvet += sym.expand((a * x + b - y) ** 2)
    return otvet

def f2(a, b):
    global dannie
    otvet = 0
    for i in range(len(dannie)):
        x = dannie[i][0]
        y = dannie[i][1]
        otvet += sym.expand(abs(a * x + b - y))
    return otvet


def dihotomia(func, a, b, delta):
    i = 0
    while b - a > eps:
        c = (a + b) / 2 - delta
        d = (a + b) / 2 + delta
        f1 = func(c)
        f2 = func(d)
        if f1 < f2:
            b = d
        else:
            a = c
        i += 1
    return (a + b) / 2




def coordinate_spusk(f, x, y):
    global eps
    worksheet = workbook.add_worksheet()
    worksheet.write(0, 0, 'n')
    worksheet.write(0, 1, 'X(n)')
    worksheet.write(0, 2, '||X(n)-X(n-1)||')
    worksheet.write(0, 3, 'f(X(n))')
    worksheet.write(0, 4, '|| f(X(n)) - f(X(n-1)) ||')
    min_x = x
    min_y = y
    del_f = 10
    i = 1
    while del_f > 2 * eps:
        del_f = f(min_x, min_y)
        tochka = [min_x, min_y]
        f_x = lambda x: f(x, min_y)
        min_x = dihotomia(f_x, -10, 10, delta)
        del_f = abs(del_f - f(min_x, min_y))
        shag = np.linalg.norm([min_x - tochka[0], min_y - tochka[1]])
        worksheet.write(i, 0, i)
        worksheet.write(i, 1, "({:.6f}, {:.6f})".format(min_x, min_y))
        worksheet.write(i, 2, "{:.6f}".format(shag))
        worksheet.write(i, 3, "{:.6f}".format(f(min_x, min_y)))
        worksheet.write(i, 4, "{:.6f}".format(del_f))
        i += 1
        if del_f < 2 * eps:
            break
        del_f = f(min_x, min_y)
        tochka = [min_x, min_y]
        f_y = lambda y: f(min_x, y)
        min_y = dihotomia(f_y, -10, 10, delta)
        del_f = abs(del_f - f(min_x, min_y))
        shag = np.linalg.norm([min_x - tochka[0], min_y - tochka[1]])
        worksheet.write(i, 0, i)
        worksheet.write(i, 1, "({:.6f}, {:.6f})".format(min_x, min_y))
        worksheet.write(i, 2, "{:.6f}".format(shag))
        worksheet.write(i, 3, "{:.6f}".format(f(min_x, min_y)))
        worksheet.write(i, 4, "{:.6f}".format(del_f))
        i += 1


coordinate_spusk(f1, x0, y0)
coordinate_spusk(f2, x0, y0)
workbook.close()
