#
#  NAME
#    problem_set1.py
#
#  DESCRIPTION
#    Open, view, and analyze raw extracellular data
#    In Problem Set 1, you will write create and test your own spike detector.
#

import numpy as np
import matplotlib.pylab as plt

def load_data(filename):
    """
    load_data takes the file name and reads in the data.  It returns two
    arrays of data, the first containing the time stamps for when they data
    were recorded (in units of seconds), and the second containing the
    corresponding voltages recorded (in units of microvolts - uV)
    """
    data = np.load(filename)[()];
    return np.array(data['time']), np.array(data['voltage'])

def bad_AP_finder(time,voltage):
    """
    This function takes the following input:
        time - vector where each element is a time in seconds
        voltage - vector where each element is a voltage at a different time

        We are assuming that the two vectors are in correspondance (meaning
        that at a given index, the time in one corresponds to the voltage in
        the other). The vectors must be the same size or the code
        won't run

    This function returns the following output:
        APTimes - all the times where a spike (action potential) was detected

    This function is bad at detecting spikes!!!
        But it's formated to get you started!
    """

    #Let's make sure the input looks at least reasonable
    if (len(voltage) != len(time)):
        print "Can't run - the vectors aren't the same length!"
        APTimes = []
        return APTimes

    numAPs = np.random.randint(0,len(time))//10000 #and this is why it's bad!!

    # Now just pick 'numAPs' random indices between 0 and len(time)
    APindices = np.random.randint(0,len(time),numAPs)

    # By indexing the time array with these indices, we select those times
    APTimes = time[APindices]

    # Sort the times
    APTimes = np.sort(APTimes)

    return APTimes

def good_AP_finder(time,voltage):
    """
    This function takes the following input:
        time - vector where each element is a time in seconds
        voltage - vector where each element is a voltage at a different time

        We are assuming that the two vectors are in correspondance (meaning
        that at a given index, the time in one corresponds to the voltage in
        the other). The vectors must be the same size or the code
        won't run

    This function returns the following output:
        APTimes - all the times where a spike (action potential) was detected
    """

    APTimes = []
    voltage = smooth(voltage, window='blackman' , window_len=9)
    #Let's make sure the input looks at least reasonable
    if (len(voltage) != len(time)):
        print "Can't run - the vectors aren't the same length!"
        return APTimes

    all_max_indices = (np.diff(np.sign(np.diff(voltage))) < 0).nonzero()[0] + 1
    all_min_indices = (np.diff(np.sign(np.diff(voltage))) > 0).nonzero()[0] + 1

    print all_min_indices, all_max_indices


    v_m = voltage.max()
    v_std = voltage.std()
    num_stds = (v_m/v_std)
    triger = (v_std * 3)
    #print "TRIGER: ", triger
    #print "MAX: ", v.max()
    #print "MEDIAN:", np.median(v)
    #print "ACG: ", np.average(v)

    triger = ((voltage.max() + triger)/2) - (voltage.max()*0.26) #this need to be improved with outliers analysis
    #print "value: ", triger
    min_indices = [time[i] for i in all_min_indices if voltage[i] < -triger]
    max_indices = [time[i] for i in all_max_indices if voltage[i] >  triger]



    #print "MAX:", max_indices
    if len(max_indices) != len(min_indices):
        max_indices = [time[i] for i in all_max_indices if voltage[i] >  triger] #118
        APTimes = max_indices
    else:
        APTimes = np.sort([(min_indices[i] + max_indices[i])/2 for i in xrange(min(len(max_indices), len(min_indices)))])


    #TODO: I can remove double detected spikes here knowing the width of the mean spike in the data
    return APTimes


def get_actual_times(dataset):
    """
    Load answers from dataset
    This function takes the following input:
        dataset - name of the dataset to get answers for

    This function returns the following output:
        APTimes - spike times
    """
    return np.load(dataset)

def detector_tester(APTimes, actualTimes):
    """
    returns percentTrueSpikes (% correct detected) and falseSpikeRate
    (extra APs per second of data)
    compares actual spikes times with detected spike times
    This only works if we give you the answers!
    """

    JITTER = 0.025 #2 ms of jitter allowed

    #first match the two sets of spike times. Anything within JITTER_MS
    #is considered a match (but only one per time frame!)

    #order the lists
    detected = np.sort(APTimes)
    actual = np.sort(actualTimes)

    #remove spikes with the same times (these are false APs)
    temp = np.append(detected, -1)
    detected = detected[plt.find(plt.diff(temp) != 0)]

    #find matching action potentials and mark as matched (trueDetects)
    trueDetects = [];
    for sp in actual:
        z = plt.find((detected >= sp-JITTER) & (detected <= sp+JITTER))
        if len(z)>0:
            for i in z:
                zz = plt.find(trueDetects == detected[i])
                if len(zz) == 0:
                    trueDetects = np.append(trueDetects, detected[i])
                    break;
    percentTrueSpikes = 100.0*len(trueDetects)/len(actualTimes)

    #everything else is a false alarm
    totalTime = (actual[len(actual)-1]-actual[0])
    falseSpikeRate = (len(APTimes) - len(actualTimes))/totalTime

    print 'Action Potential Detector Performance performance: '
    print '     Correct number of action potentials = ' + str(len(actualTimes))
    print '     Percent True Spikes = ' + str(percentTrueSpikes)
    print '     False Spike Rate = ' + str(falseSpikeRate) + ' spikes/s'
    print
    return {'Percent True Spikes':percentTrueSpikes, 'False Spike Rate':falseSpikeRate}


def plot_spikes(time,voltage,APTimes,titlestr, actualTimes):
    """
    plot_spikes takes four arguments - the recording time array, the voltage
    array, the time of the detected action potentials, and the title of your
    plot.  The function creates a labeled plot showing the raw voltage signal
    and indicating the location of detected spikes with red tick marks (|)
    """

    plt.figure()

    plt.xlabel("Time (s)")
    plt.ylabel("Voltage (uV)")
    plt.title(titlestr)

    plt.plot(time,voltage)
    for apt in APTimes:
        drop_pin_at(apt)

    #for act in actualTimes:
    #    drop_pin_at(act, color='g', ymin=0.98, ymax=1.0)
    plt.show()

def plot_waveforms(time,voltage,APTimes,titlestr):
    """
    plot_waveforms takes four arguments - the recording time array, the voltage
    array, the time of the detected action potentials, and the title of your
    plot.  The function creates a labeled plot showing the waveforms for each
    detected action potential
    """

    plt.figure()
    plt.xlabel("Time (s)")
    plt.ylabel("Voltage (uV)")
    plt.title(titlestr)

    delta = np.float32(0.003)

    for apt in APTime:
        data = v[(t>apt-delta) & (t<apt+delta)]
        _time = np.linspace(-0.003, 0.003, data.shape[0])
        plt.plot(_time, data, color='b')

    plt.show()



def drop_pin_at(x, linewidth=1, color='r', ymin=0.95, ymax=0.97):
    """
    Drop a vertical red pin at the given x cord.
    Here we are dropping pins at each AP.
    """
    plt.axvline(x=x, linewidth=linewidth, color=color, ymin=ymin, ymax=ymax)



def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y[(window_len/2):-(window_len/2)]


##########################
#You can put the code that calls the above functions down here
if __name__ == "__main__":
    print "m..."
    t,v = load_data('spikes_hard_test.npy')

    actualTimes = get_actual_times('spikes_hard_practice_answers.npy')
    APTime = good_AP_finder(t,v)
    plot_spikes(t,v,APTime,'Action Poteltials in Raw Signal', actualTimes)
    plot_waveforms(t,v,APTime,'Waveforms')
    detector_tester(APTime,actualTimes)
