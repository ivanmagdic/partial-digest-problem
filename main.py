import itertools
import operator
import sys
import math
import time
from functools import reduce
import multiprocessing
from functools import partial


# x - the set of integers representing the location of all cuts in the restriction map, including the start and end L

# L - The multiset of integers representing lengths of each of the DNA fragments produced from a partial digest; formed
# from X by taking all pairwise differences.


def find_all(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1:
            return
        yield start
        start += 1


def quadratic_equation(i):
    a = -1
    b = 1
    c = 2 * i

    x = int(round(b ** 2) - (4 * a * c))

    solution_1 = (-b + math.sqrt(x)) / (2 * a)
    solution_2 = (-b - math.sqrt(x)) / (2 * a)

    if solution_1 > 0:
        return int(solution_1)

    return int(solution_2)


def count_iterable(i):
    return sum(1 for e in i)


def create_multiset2(x):
    aux_x = []
    for x in list(itertools.combinations(x, 2)):
        aux_s = reduce(operator.__sub__, x) * -1
        aux_x.append(aux_s)
    return aux_x


# Used for large lists (where n > 7 )
def another_brute_force_pdp(L, n):
    global last
    last = max(L)

    pool = multiprocessing.Pool()

    chunk = 6

    c = count_iterable(itertools.combinations(L, n - 2))

    if c < chunk:
        chunk = c

    start = 0
    stop = chunk

    while c > 0:
        pool.map(partial(process_data, L), itertools.islice(itertools.combinations(L, n - 2), start, stop))

        if c < chunk:
            chunk = c
            continue

        c -= chunk
        start += chunk
        stop += chunk

    pool.close()
    pool.join()


def brute_force_pdp(L, n):
    global last
    last = max(L)

    pool = multiprocessing.Pool(multiprocessing.cpu_count())
    pool.map(process_data, list(itertools.combinations(L, n - 2)))
    pool.close()
    pool.join()


def process_data(L, comb):
    aux_x = [0]
    aux_x += comb
    aux_x.append(last)
    # aux_x.sort()
    aux_multiset = create_multiset2(aux_x)
    aux_multiset.sort()
    if aux_multiset == L:
        print('Found: ', aux_x)
        # return aux_x


def partial_digest(L):
    global width
    width = max(L)
    L.remove(width)
    x = [0, width]
    place(L, x)


def place(L, x):
    if not L:
        aux_x = [j for j in x]
        aux_x.sort()
        print(aux_x)
        return

    y = max(L)

    if is_subset(y, x, L):
        x.append(y)
        remove_elements(y, x, L)
        place(L, x)
        if y in x:
            x.remove(y)
        L.extend(delete(y, x))

    if is_subset(abs(width - y), x, L):
        x.append(abs(width - y))
        remove_elements(abs(width - y), x, L)
        place(L, x)
        if abs(width - y) in x:
            x.remove(abs(width - y))
        L.extend(delete(abs(width - y), x))

    return


def delete(y, X):
    aux_l = []
    for i in X:
        aux_l.append(abs(y - i))
    return aux_l


def remove_elements(y, X, L):
    for xi in X:
        if abs(y - xi) in L:
            L.remove(abs(y - xi))


def is_subset(y, X, L):
    for xi in X:
        if abs(y - xi) not in L:
            return False
    return True


if __name__ == "__main__":

    dnk = open(sys.argv[1]).read()
    search = sys.argv[2].split(',')

    x = [0]

    for s in search:
        x += find_all(dnk, s)

    x = list(set(x))
    x.sort()
    x.append(len(dnk) - 1)

    L = create_multiset2(x)
    L.sort()

    n = quadratic_equation(len(L))

    print('Search value:', search)
    print('SET (množica):', n - 2)
    print('\n')

    print('MULTISET (multimnožica):')
    print(L)
    print('\n')

    diffs = []
    repetitions = 1

    # TEST
    # L = [2, 2, 3, 3, 4, 5, 6, 7, 8, 10]

    for i in range(0, repetitions):
        print('Brute force ', i)
        print('Solutions:')

        start = time.time()
        another_brute_force_pdp(L, int(n))
        # brute_force_pdp(L, int(n))
        end = time.time()

        diff = (end - start) * 1000
        diffs.append(diff)

        print('Elapsed time:', diff, 'ms')

        print('\n')

    print('Average time for', repetitions, 'repetitions is', sum(diffs) / len(diffs), 'ms')

    for i in range(0, repetitions):
        print('Partial digest', i)
        print('Solutions:')

        start = time.time()
        partial_digest(L)
        end = time.time()

        diff = (end - start) * 1000
        diffs.append(diff)

        print('Elapsed time:', diff, 'ms')

        print('\n')

    print('Average time for', repetitions, 'repetitions is', sum(diffs) / len(diffs), 'ms')
