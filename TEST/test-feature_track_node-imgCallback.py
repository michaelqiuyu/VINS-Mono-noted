


def test(frep, ffrep):
    i = 0
    k = 1
    per_count = 0
    last_count = 0

    while True:
        if i > 100000:
            return False
        per_count += 1
        if i == 0:
            i += 1
            continue

        temp = (k * frep) / i
        if ((k * frep) / i) < ffrep + 0.5:
            print("i = " + str(i) + ", k = " + str(k) + "count = " + str(per_count) + "temp = " + str(temp) + "frep = " + str(frep) + "ffrep = " + str(ffrep))
            if last_count == 0:
                last_count = per_count

            if abs(temp - ffrep) < 0.01 * ffrep:
                if i > 100 and abs(last_count - per_count) > 5:
                    print("找到i = " + str(i) + ", k = " + str(k) + ", frep = " + str(frep) + "ffrep = " + str(ffrep) + "last = " + str(last_count) + "per_count = " + str(per_count))
                    return True

            last_count = per_count
            per_count = 0

            k += 1
            i += 1
        else:
            i += 1


def testt():
    for ffrep in range(5, 30):
        for frep in range(1 + ffrep, ffrep + 100):
            if test(frep, ffrep):
                return


if __name__ == '__main__':
    testt()



