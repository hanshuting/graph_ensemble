import pylab

def create_array_from_lines(lines):
    l = len(lines)
    array = [0]*l
    i = 0
    for line in lines:
        array[i] = eval(line)
        i += 1
    return array


xilines        = open("xi.dump").readlines()
losslines      = open("loss.dump").readlines()
regvallines    = open("regval.dump").readlines()
exactobjlines  = open("exactobj.dump").readlines()
approxobjlines = open("approxobj.dump").readlines()
epsilonlines   = open("epsilon.dump").readlines()
gammalines     = open("gamma.dump").readlines()

epsilon   = create_array_from_lines(epsilonlines)
xi        = create_array_from_lines(xilines)
gamma     = create_array_from_lines(gammalines)
regval    = create_array_from_lines(regvallines)
exactobj  = create_array_from_lines(exactobjlines)
approxobj = create_array_from_lines(approxobjlines)
loss      = create_array_from_lines(losslines)

pylab.figure()

pylab.subplot(311)
pylab.title("loss")
pylab.plot(loss)

pylab.subplot(312)
pylab.title("xi")
pylab.plot(xi)

pylab.subplot(313)
pylab.title("regularizer value")
pylab.plot(regval)

pylab.figure()
pylab.subplot(311)
pylab.title("exact obj val")
pylab.semilogy(exactobj)

pylab.subplot(312)
pylab.title("approx obj val")
pylab.semilogy(approxobj)

pylab.subplot(313)
pylab.title("epsilon: exact obj val - approx obj val")
pylab.semilogy(epsilon)

pylab.figure()
pylab.title("gamma: min. exact obj bal - approx obj val")
pylab.semilogy(gamma)

pylab.show()


