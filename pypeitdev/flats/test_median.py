
import numpy as np
from matplotlib import pyplot as plt
# Imports for this fast running median routine
from collections import deque
from itertools import islice
from bisect import insort, bisect_left



import scipy.ndimage

def running_median_scipy_medfilt(seq, window_size):
    return scipy.ndimage.median_filter(seq, size = window_size,mode='reflect')
#    return (seq, window_size)
    #if window_size % 2 == 1 else window_size + 1)
#    return medfilt(seq, window_size if window_size % 2 == 1 else window_size + 1)
#    return medfilt(seq, window_size)


def running_median_insort(seq, window_size):
    """Contributed by Peter Otten"""

    # pad the array for the reflection
    seq_pad = np.concatenate((seq[0:window_size][::-1],seq,seq[-1:(-1-window_size):-1]))

    window_size= int(window_size)
    seq_pad = iter(seq_pad)
    d = deque()
    s = []
    result = []
    for item in islice(seq_pad, window_size):
        d.append(item)
        insort(s, item)
        result.append(s[len(d)//2])
    m = window_size // 2
    for item in seq_pad:
        old = d.popleft()
        d.append(item)
        del s[bisect_left(s, old)]
        insort(s, item)
        result.append(s[m])

    result = np.roll(result, -window_size//2 + 1)
    return result[window_size:-window_size]


window = 1001
sequence = np.random.rand(50*window)


out_scipy = running_median_scipy_medfilt(sequence, window)
out_insort= running_median_insort(sequence,window)

plt.plot(out_scipy - out_insort)
#plt.plot(out_insort)
#plt.plot(np.roll(out_insort,-window//2 + 1))
plt.show()