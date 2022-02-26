from Initial_Values import *
import math
import matplotlib.pyplot as plt
import numpy as np
# TIME_DOMAIN
def time_vector(NUM_TOTAL_SAMPLES,SAMPLE_TIME_S):
    timeVector = []
    for x in range(NUM_TOTAL_SAMPLES):
        y = x * SAMPLE_TIME_S
        timeVector.append(y)
    return timeVector
def ComputeResponse(NUM_TOTAL_SAMPLES,CUTOFF_FREQUENCY_HZ,CUTOFF_FREQUENCY2_HZ,SAMPLE_TIME_S,NUM_SHIFT_SAMPLES,choice):
    impulseResponse = []
    PI = math.pi
    for n in range(NUM_TOTAL_SAMPLES):
        if n != NUM_SHIFT_SAMPLES:
            if choice == 1:
                y = math.sin(2.0 * PI * CUTOFF_FREQUENCY_HZ * SAMPLE_TIME_S * (n - NUM_SHIFT_SAMPLES)) / (PI * SAMPLE_TIME_S * (n - NUM_SHIFT_SAMPLES))
                impulseResponse.append(y)
            if choice == 2:
                y = (math.sin(PI * (n - NUM_SHIFT_SAMPLES)) - math.sin(2.0 * PI * CUTOFF_FREQUENCY_HZ * SAMPLE_TIME_S * (n - NUM_SHIFT_SAMPLES))) / (PI * SAMPLE_TIME_S * (n - NUM_SHIFT_SAMPLES))
                impulseResponse.append(y)
            if choice == 3:
                y = (math.sin(2.0 * PI * CUTOFF_FREQUENCY2_HZ * SAMPLE_TIME_S * (n - NUM_SHIFT_SAMPLES)) - math.sin(2.0 * PI * CUTOFF_FREQUENCY_HZ * SAMPLE_TIME_S * (n - NUM_SHIFT_SAMPLES))) / (PI * SAMPLE_TIME_S * (n - NUM_SHIFT_SAMPLES))
                impulseResponse.append(y)
            if choice == 4:
                y = (math.sin(2.0 * PI * CUTOFF_FREQUENCY_HZ * SAMPLE_TIME_S * (n - NUM_SHIFT_SAMPLES)) - math.sin(2.0 * PI * CUTOFF_FREQUENCY2_HZ * SAMPLE_TIME_S * (n - NUM_SHIFT_SAMPLES)) + math.sin(PI * (n - NUM_SHIFT_SAMPLES))) / (PI * SAMPLE_TIME_S * (n - NUM_SHIFT_SAMPLES))
                impulseResponse.append(y)
        else:
            if choice == 1:
                y = 2.0 * CUTOFF_FREQUENCY_HZ
                impulseResponse.append(y)
            if choice == 2:
                y = 1.0 / SAMPLE_TIME_S - 2.0 * CUTOFF_FREQUENCY_HZ
                impulseResponse.append(y)
            if choice == 3:
                y = 2.0 * CUTOFF_FREQUENCY2_HZ - 2.0 * CUTOFF_FREQUENCY_HZ
                impulseResponse.append(y)
            if choice == 4:
                y = 2.0 * CUTOFF_FREQUENCY_HZ - 2.0 * CUTOFF_FREQUENCY2_HZ + 1.0 / SAMPLE_TIME_S
                impulseResponse.append(y)
    for n in range(NUM_TOTAL_SAMPLES):
        impulseResponse[n] *= SAMPLE_TIME_S
    return impulseResponse
def ComputeWindow(NUM_TOTAL_SAMPLES,choice2):
    window= []
    PI = math.pi
    for n in range(NUM_TOTAL_SAMPLES):
        if choice2 == 1:
            y = 1
            window.append(y)
        if choice2 == 2:
            y = 1.0 - abs((n - 0.5 * NUM_TOTAL_SAMPLES) / (0.5 * NUM_TOTAL_SAMPLES))
            window.append(y)
        if choice2 == 3:
            y = 1.0 - pow((n - 0.5 * NUM_TOTAL_SAMPLES) / (0.5 * NUM_TOTAL_SAMPLES), 2.0)
            window.append(y)
        if choice2 == 4:
            y = math.sin(PI * n / NUM_TOTAL_SAMPLES)
            window.append(y)
        if choice2 == 5:
            y = 0.5 * (1 - math.cos(2.0 * PI * n / NUM_TOTAL_SAMPLES))
            window.append(y)
        if choice2 == 6:
            y = (25.0 / 46.0) - (21.0 / 46.0) * math.cos(2.0 * PI * n / NUM_TOTAL_SAMPLES)
            window.append(y)
        if choice2 == 7:
            y = 0.42 - 0.5 * math.cos(2.0 * PI * n / NUM_TOTAL_SAMPLES) + 0.08 * math.cos(4.0 * PI * n / NUM_TOTAL_SAMPLES)
            window.append(y)
        if choice2 == 8:
            y = 0.355768 - 0.487396 * math.cos(2.0 * PI * n / NUM_TOTAL_SAMPLES) + 0.144232 * math.cos(4.0 * PI * n / NUM_TOTAL_SAMPLES) - 0.012604 * math.cos(6.0 * PI * n / NUM_TOTAL_SAMPLES)
            window.append(y)
        if choice2 == 9:
            y = 0.3635819 - 0.4891775 * math.cos(2.0 * PI * n / NUM_TOTAL_SAMPLES) + 0.1365995 * math.cos(4.0 * PI * n / NUM_TOTAL_SAMPLES) - 0.0106411 * math.cos(6.0 * PI * n / NUM_TOTAL_SAMPLES)
            window.append(y)
        if choice2 == 10:
            y = 0.35875 - 0.48829 * math.cos(2.0 * PI * n / NUM_TOTAL_SAMPLES) + 0.14128 * math.cos(4.0 * PI * n / NUM_TOTAL_SAMPLES) - 0.01168 * math.cos(6.0 * PI * n / NUM_TOTAL_SAMPLES)
            window.append(y)
        if choice2 == 11:
            y = 0.21557895 - 0.41663158 * math.cos(2.0 * PI * n / NUM_TOTAL_SAMPLES) + 0.277263158 * math.cos(4.0 * PI * n / NUM_TOTAL_SAMPLES) - 0.083578947 * math.cos(6.0 * PI * n / NUM_TOTAL_SAMPLES) + 0.006947368 * math.cos(8.0 * PI * n / NUM_TOTAL_SAMPLES)
            window.append(y)
    return window
def WindowedResponse(NUM_TOTAL_SAMPLES,impulseResponse,ComputeWindow):
    windowedImpulseResponse = []
    for n in range(NUM_TOTAL_SAMPLES):
        y = impulseResponse[n] * ComputeWindow[n]
        windowedImpulseResponse.append(y)
    return windowedImpulseResponse
#FREQ_DOMAIN
def Freq_vector(NUM_FREQ_SAMPLES):
    frequencyvector = []
    df = (0.5 / SAMPLE_TIME_S) / (NUM_FREQ_SAMPLES - 1.0)
    for n in range(NUM_FREQ_SAMPLES):
        y = n * df
        frequencyvector.append(y)
    return frequencyvector
def ComputeRespBode_Window(NUM_FREQ_SAMPLES,NUM_TOTAL_SAMPLES,Freq_vector,SAMPLE_TIME_S,WindowedResponse):
    winRespMag = []
    for n in range(NUM_FREQ_SAMPLES):
        reWin = 0.0
        imWin = 0.0
        for x in range(NUM_TOTAL_SAMPLES):
            reWin = reWin + WindowedResponse[x] * math.cos(2.0 * math.pi * Freq_vector[n] * x * SAMPLE_TIME_S)
            imWin = imWin + WindowedResponse[x] * math.sin(2.0 * math.pi * Freq_vector[n] * x * SAMPLE_TIME_S)
        d = 10.0 * math.log10(reWin * reWin + imWin * imWin)
        winRespMag.append(d)
    return winRespMag
def ComputeRespBode_NoWindow(NUM_FREQ_SAMPLES,NUM_TOTAL_SAMPLES,ComputeResponse,Freq_vector,SAMPLE_TIME_S):
    impRespMag = []
    for n in range(NUM_FREQ_SAMPLES):
        re = 0.0
        im = 0.0
        for x in range(NUM_TOTAL_SAMPLES):
            re = re + ComputeResponse[x] * math.cos(2.0 * math.pi * Freq_vector[n] * x * SAMPLE_TIME_S)
            im = im + ComputeResponse[x] * math.sin(2.0 * math.pi * Freq_vector[n] * x * SAMPLE_TIME_S)
        y = 10.0 * math.log10(re * re + im * im)
        impRespMag.append(y)
    return impRespMag
def ComputeWinDFT(NUM_FREQ_SAMPLES,NUM_TOTAL_SAMPLES,ComputeWindow,Freq_vector,SAMPLE_TIME_S):
    WinDTF = []
    for x in range(NUM_FREQ_SAMPLES):
        re = 0.0
        im = 0.0
        for n in range(NUM_TOTAL_SAMPLES):
            re = re + ComputeWindow[n] * math.cos(2.0 * math.pi * Freq_vector[x] * n * SAMPLE_TIME_S)
            im = im + ComputeWindow[n] * math.sin(2.0 * math.pi * Freq_vector[x] * n * SAMPLE_TIME_S)
        y = 10 * math.log10(re * re + im * im)
        WinDTF.append(y)
    return WinDTF

a = ComputeResponse(NUM_TOTAL_SAMPLES,CUTOFF_FREQUENCY_HZ,CUTOFF_FREQUENCY2_HZ,SAMPLE_TIME_S,NUM_SHIFT_SAMPLES,Filter_type)
i = time_vector(NUM_TOTAL_SAMPLES,SAMPLE_TIME_S)
v = Freq_vector(NUM_FREQ_SAMPLES)
f = ComputeWindow(NUM_TOTAL_SAMPLES,Window_type)
d = WindowedResponse(NUM_TOTAL_SAMPLES,a,f)
x = ComputeRespBode_Window(NUM_FREQ_SAMPLES,NUM_TOTAL_SAMPLES,v,SAMPLE_TIME_S,d)
n = ComputeRespBode_NoWindow(NUM_FREQ_SAMPLES,NUM_TOTAL_SAMPLES,a,v,SAMPLE_TIME_S)
r = ComputeWinDFT(NUM_FREQ_SAMPLES,NUM_TOTAL_SAMPLES,f,v,SAMPLE_TIME_S)



fig, axs = plt.subplots(2, 2)
axs[0, 0].plot(i,a)
axs[0, 0].set_title('Impulse Response')
axs[0, 1].plot(i, f, 'tab:orange')
axs[0, 1].set_title('Window response')
axs[1, 0].plot(v, x, 'tab:green')
axs[1, 0].plot(v, n, 'tab:red')
axs[1, 0].set_title('Windowed vs unwindowed DTFs')
axs[1, 1].plot(v, r, 'tab:red')
axs[1, 1].set_title('Window DFT')
plt.show()