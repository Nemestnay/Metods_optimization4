import numpy as np
import sympy as sym
import xlsxwriter
from scipy.misc import derivative


dannie = [[1, 0], [2, -1], [3, -3], [4, 2], [5, -2]]
eps = 0.000001
delta = eps / 4
workbook = xlsxwriter.Workbook('КР4_2.xlsx')
x0 = 10
y0 = 10
a_min, a_max, b_min, b_max = -1, 0.5, -1, 1
L = 70

def f1(a, b):
    return sum([pow(a * dannie[k][0] + b - dannie[k][1], 2) for k in range(5)])

def grad(f, a, b):
    return [derivative(lambda x: f(x, b), a, dx=eps), derivative(lambda y: f(a, y), b, dx=eps)]

def grad_f1(a, b):
    return grad(f1, a, b)


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

def vibor_shaga(x, y, f, grad):
    delta = 0.25
    alpha = 1
    while f(x-alpha*grad(x,y)[0], y-alpha*grad(x,y)[1]) > f(x,y) - delta*alpha*np.linalg.norm(grad(x, y))**2:
        alpha /= 2
    return alpha


def gm_droblenie_shaga(f, grad, x0, y0):
    global eps
    worksheet = workbook.add_worksheet()
    worksheet.write(0, 0, 'градиентный метод с дроблением шага')
    worksheet.write(1, 0, 'n')
    worksheet.write(1, 1, 'X(n)')
    worksheet.write(1, 2, 'f(X(n))')
    worksheet.write(1, 3, "|| f'(X(n)) ||")
    min_x = x0
    min_y = y0
    del_f = 10
    i = 0
    while del_f > eps and i < 1000:
        step = vibor_shaga(min_x, min_y, f, grad)
        eps1 = f(min_x, min_y)
        dx = -step * grad(min_x, min_y)[0]
        dy = -step * grad(min_x, min_y)[1]
        min_x += dx
        min_y += dy
        del_f = abs(eps1 - f(min_x, min_y))
        shag = np.linalg.norm(grad(min_x, min_y))
        if i % 10 == 0:
            worksheet.write(i//10 +2, 0, i)
            worksheet.write(i//10+2, 1, "({:.6f}, {:.6f})".format(min_x, min_y))
            worksheet.write(i//10+2, 2, "{:.6f}".format(f(min_x, min_y)))
            worksheet.write(i//10+2, 3, "{:.6f}".format(shag))
        i += 1


def gm_postoyanniy_shag(f, grad, x0, y0, step):
    global eps
    worksheet = workbook.add_worksheet()
    worksheet.write(0, 0, 'градиентный метод с постоянным шагом')
    worksheet.write(1, 0, 'n')
    worksheet.write(1, 1, 'X(n)')
    worksheet.write(1, 2, 'f(X(n))')
    worksheet.write(1, 3, "|| f'(X(n)) ||")
    min_x = x0
    min_y = y0
    i = 1
    del_f = 10
    while del_f > eps and i < 1000:
        del_f = f(min_x, min_y)
        dx = -step * grad(min_x, min_y)[0]
        dy = -step * grad(min_x, min_y)[1]
        min_x += dx
        min_y += dy
        del_f = abs(del_f - f(min_x, min_y))
        shag = np.linalg.norm(grad(min_x, min_y))
        if i % 10 == 0:
            worksheet.write(i//10 +1, 0, i)
            worksheet.write(i//10+1, 1, "({:.6f}, {:.6f})".format(min_x, min_y))
            worksheet.write(i//10+1, 2, "{:.6f}".format(f(min_x, min_y)))
            worksheet.write(i//10+1, 3, "{:.6f}".format(shag))
        i += 1

def gm_naiscoreyshiy_spusk(f, grad, x0, y0):
    global eps
    worksheet = workbook.add_worksheet()
    worksheet.write(0, 0, 'наискорейший градиентный спуск')
    worksheet.write(1, 0, 'n')
    worksheet.write(1, 1, 'X(n)')
    worksheet.write(1, 2, 'f(X(n))')
    worksheet.write(1, 3, "|| f'(X(n)) ||")
    min_x = x0
    min_y = y0
    del_f = 10
    i = 0
    while del_f > eps and i < 1000:
        step = dihotomia(lambda x: f(min_x-x*grad(min_x,min_y)[0], min_y-x*grad(min_x,min_y)[1]), 0, 100, eps / 10)
        del_f = f(min_x, min_y)
        dx = -step * grad(min_x, min_y)[0]
        dy = -step * grad(min_x, min_y)[1]
        min_x += dx
        min_y += dy
        del_f = abs(del_f - f(min_x, min_y))
        shag = np.linalg.norm(grad(min_x, min_y))
        worksheet.write(i+2, 0, i)
        worksheet.write(i+2, 1, "({:.6f}, {:.6f})".format(min_x, min_y))
        worksheet.write(i+2, 2, "{:.6f}".format(f(min_x, min_y)))
        worksheet.write(i+2, 3, "{:.6f}".format(shag))
        i += 1


gm_droblenie_shaga(f1, grad_f1, x0, y0)
gm_postoyanniy_shag(f1, grad_f1, x0, y0, 1 / L)
gm_naiscoreyshiy_spusk(f1, grad_f1, x0, y0)
workbook.close()
