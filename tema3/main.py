def citire(filename):
    with open(filename, 'r') as f:
        n = int(f.readline())
        b = []
        a = [[] for _ in range(n)]
        for line in f.readlines():
            if line == '\n':
                continue
            val = float(line.split(',')[0])  # luam valoarea
            if len(line.split(',')) > 1:
                i = int(line.split(',')[1])
                j = int(line.split(',')[2])
                nuExista = True
                for l in range(len(a[i])):
                    # mai multe valori cu aceiasi indici de linie si coloana
                    if a[i][l][1] == j:
                        a[i][l][0] += val
                        nuExista = False
                        break
                if nuExista:
                    a[i].append([val, j])
            else:
                b.append(val)
    return a, b


def egale(a, b):
    esp = 1e-16
    for i in range(len(a)):
        for val, j in b[i]:
            for l in range(len(a[i])):
                if a[i][l][1] == j:
                    dif = a[i][l][0] - val
                    if dif > esp:
                        return False
    return True


def aduna(a, b):
    n = len(a)
    rez = a
    for i in range(n):
        for val, j in b[i]:
            nuExista = True
            for l in range(len(rez[i])):
                if rez[i][l][1] == j:
                    rez[i][l][0] += val
                    nuExista = False
                    break
            if nuExista:
                rez[i].append([val, j])
    return rez


def inmultesteMatrici(a, b):
    mat = []
    for i in range(len(a)):
        linie = []
        for j in range(len(b)):
            sum = 0
            for l in range(len(a[i])):
                bb = 0
                for x in range(len(b[a[i][l][1]])):
                    if b[a[i][l][1]][x][1] == j % (len(b)):
                        bb = b[a[i][l][1]][x][0]
                        break
                sum = sum + bb * a[i][l][0]
            if sum != 0:
                linie.append([sum, j])
        mat.append(linie)
    return mat


def start():
    # adunare
    a = citire('a.txt')[0]
    b = citire('b.txt')[0]
    aplusb = citire('a_plus_b.txt')[0]
    a = aduna(a, b)
    isEqual1 = egale(a, aplusb)
    print("A + B = a_plus_b: ", isEqual1)

    # inm
    a = citire('a.txt')[0]
    aoria = citire('a_ori_a.txt')[0]
    a = inmultesteMatrici(a, a)
    isEqual2 = egale(a, aoria)
    print("A * A = a_ori_a:", isEqual2)


if __name__ == '__main__':
    start()
