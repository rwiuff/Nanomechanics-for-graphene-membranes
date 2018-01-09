import matplotlib.pyplot as plt
import numpy as np
fig = plt.figure()
x = np.zeros(5)
y = np.zeros(5)
RMS = np.array([[0.34876869,  0.18938546, 0.12407506, 0.10074879,  0.08974112],
                [0.81348904,  0.42539935, 0.25613696, 0.21323875,  0.18362789],
                [1.62905479,  0.82789468, 0.49859411, 0.42206707,  0.35481848],
                [2.72562078,  1.36543878, 0.82432755, 0.70459229,  0.5859104],
                [3.86834624,  1.86860316, 1.14843955, 1.0077887,  0.82560262],
                [5.76951934,  2.84784375, 1.7247009, 1.48853344,  1.22645908],
                [7.7216739,  3.79402227, 2.29978961, 1.98998895,  1.63601802],
                [9.96458781,  4.87791536, 2.95855463, 2.56468172,  2.1052925],
                [12.50183139, 6.10047193, 3.70136748, 3.21287784,   2.63447271],
                [15.33748917,   7.46276768, 4.52865011, 3.93489,   3.22377782]])
Index = [0, 1, 3, 5, 6]
for i in range(10):
    plt.plot(Index, RMS[i, :], label='{} nm'.format(i + 1))
    plt.legend(loc=1)
    plt.ylabel('RMS')
    kticks = Index
    ticklabels = ['0', '1', '3', '5', '6']
    plt.xticks(kticks, ticklabels)
    plt.xlabel('Modes')
    plt.title('RMS for varying single layer membrane hole sizes')
plt.subplots_adjust(left=0.05, bottom=None, right=0.98, top=0.94,
                    wspace=None, hspace=None)
plt.show()
