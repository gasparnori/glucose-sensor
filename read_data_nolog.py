import numpy as np
import matplotlib.pyplot as plt
import serial
import traceback
import binascii
from datetime import datetime, timedelta
import serial.tools.list_ports
import time
import multiprocessing as mp
import csv

# if False, it uses the settings of the IVS computer
Debugging = False
TIMEthres=1000

if not Debugging:
    ###################### defining paths on the experiment computer ###########
    saveOutput = True
    filePath = 'C:/GluSense/Reference measurements/'  # where to save the output
    Measurement_Rate = 1  ## in seconds
    calibDataPath = 'C:/GluSense/jobst_glucose_sensor/device calibration data/259 1R.def'
    logPath='C:/GluSense/GluSense Monitoring Software/Log/InVitroApp.txt'
    # depending on the computer it can be a fixed port or not
    comPort = 'COM7'  ## needs to be adjusted for the computer
    fixedPort = False
else:
    ###################### defining paths for debugging ###########
    saveOutput = False
    filePath = ''  # where to save the output
    Measurement_Rate = 10
    logPath = 'InVitroApp.txt'
    calibDataPath = '259 1R.def'
    # depending on the computer it can be a fixed port or not
    comPort = 'COM7'
    fixedPort = False

### fields for the csv file
outputfields = ['time', 'glucose [mM]', 'glucose [mM] before temperature compensation', 'temperature [C]',
                    'temperature constant']

## low pass filter applied on the glucose data
lpf=0.9
#plotting window
plotwin=500
medianfiltwin=10

class calibrationData:
    """

    class that loads the .def file received with the device, and interprets the values in it
    It uses two possible current ranges: current_range = {'25nA': 1, '50nA': 2}

    """
    def __init__(self):
        self.currentRange=2 ##
        calibline, config, subtract, tempsens=None, None, None, None
        with open(calibDataPath, "r", encoding="utf8") as calibfile:
            [calibline, _, _, config, subtract, _, tempsens, _] = calibfile.readlines()
        calibline= [float(c) for c in calibline.split('\n')[0].split(',')]
        subtract = [int(c) for c in subtract.split('\n')[0].split(',')]
        tempsens = [float(c) for c in tempsens.split('\n')[0].split(',')]

        #the location of the relevant data in the table
        self.glucosech  = config.split(',').index('"Glucose"')

        ## the channel to subtract the glucose channel from
        #  -1 is needed because: the numbering in the configuration starts from 1,
        # the numbering in python starts from 0...
        self.subtractch = subtract[self.glucosech]-1

        ## the constant for converting the results to mM
        self.convconstant = calibline[self.glucosech]/100 ##* self.currentRange
        self.tempconstant = tempsens[self.glucosech]/100
        self.defaulttemp = tempsens[7]

    def calc_glucose(self, m):
        """

        :param m: measurement object with channel values
        :return: calculated glucose, calculated glucose with temperature compensation

        """
        gluc =((m.channels[self.glucosech]- m.channels[self.subtractch]) * self.convconstant)
        gluc_comp = gluc/(np.exp(self.tempconstant * (m.temperature - self.defaulttemp)))
        return gluc, gluc_comp


class MeasurementPoint:
    """
    An object that reads and interprets all the bytes received in one measurement
    """
    def __init__(self, data, timestamp):
        self.timestamp = timestamp
        self.timestring = datetime.strftime(self.timestamp, '%d/%m/%Y %H:%M')
        self.data = str(binascii.hexlify(data), 'ascii')
        self.data = [self.data[i:i + 2] for i in range(0, len(self.data), 2)]
        # checking if we got real data
        validity = self.check_validity()  # checks start bytes, length bytes and stop byte
        messagetype = int(self.data[4], 16)
        if messagetype == 4 and validity:
            self.validData = True
            ## reading byte number 5 to 17
            self.channels=[self.get_channel_data(i*2+5) for i in range(0,6)]
            self.temperature = (self.get_channel_data(17) / 16)
            self.id = self.get_id(19)  # from index 38 to 46
            # print (self.id)
        elif messagetype == 5 and validity:
            print('error code:', self.data[10:12])
        else:
            print('invalid data received')

    def check_validity(self):
        """

        :return: a boolean whether the start bytes, stop bytes and the checksum are correct
        """
        startbyte1 = self.data[0]
        startbyte2 = self.data[3]
        if startbyte1 != '68' or startbyte2 != '68':
            print('start byte error')
            return False
        lenbyte1 = int(self.data[1], 16)
        lenbyte2 = int(self.data[2], 16)
        if lenbyte1 != 19 or lenbyte2 != 19:
            print('len byte error')
            return False
        stopbtye = self.data[24]
        if stopbtye != '16':
            print('stop byte error')
            return False
        checksum = str(hex(np.sum([int(self.data[i], 16) for i in range(4, 23)])))[-2:]
        if checksum != self.data[23]:
            print('checksum not matching')
            return False
        return True

    def get_channel_data(self, startind):
        """

        :param startind: start index of the two bytes of channel data
        :return: data in int format
        """
        ### to do: deal with negative data if data>2^15
        ## check for errors: if the values are exactly 32767 or -32768
        highbyte = np.left_shift(int(self.data[startind], 16), 8)
        lowbyte = int(self.data[startind + 1], 16)
        data = highbyte + lowbyte
        if data>32767:  #if the data is above 2^15: it represents negative data
            data=data-65536  ## not sure if this is needed for the
        # elif data==32767 or data==-32768:
        #    print('measurement is out of range')
        if data == 65535:
            print('measurement is out of range')

        return data

    def get_id(self, startindex):
        """

        :param startindex: start index of the four bytes data
        :return: ID number
        """
        b0 = np.left_shift(int(self.data[startindex], 16), 24)
        b1 = np.left_shift(int(self.data[startindex + 1], 16), 16)
        b2 = np.left_shift(int(self.data[startindex + 2], 16), 8)
        b3 = int(self.data[startindex + 3], 16)
        return b0 + b1 + b2 + b3


class Plotter:
    """
    Plotter object to plot measurements in real time
    """
    def __init__(self):
        self.gluc=[]
        self.temp=[]
        self.times=[]
        self.gluc_smooth=[]
        self.refconc=[]
        self.f,[self.axis, self.tempaxis]=plt.subplots(2, sharex=True)
        # self.tempaxis=self.axis.twinx()
        self.axis.set_title('glucose and temperature measurements')
        self.axis.set_ylabel('glucose [mM]')
        self.tempaxis.set_ylabel('temperature [C]')
        self.tempaxis.set_xlabel('time')
        self.refaxis=self.axis.twinx()
        self.refaxis.set_ylim([0, 400])
        self.refaxis.set_ylabel('theoretical glucose [mg/dL]')
        self.templine, = self.tempaxis.plot([], [], '.b', label='temperature')
        self.glucline, = self.axis.plot([], [], '*r', alpha=0.5, label='glucose')
        self.glucline_smooth, = self.axis.plot([], [], '--', color='k', label='filtered')
        self.theoretical_line, =self.refaxis.plot([],[], 'g', label='theoretical glucose')
        self.f.legend()

    def update_plot(self,time, gluc, temp, gluc_smooth, theoretical_conc):
        """

        :param time:
        :param gluc:
        :param temp:
        :param gluc_smooth:
        :param theoretical_conc:
        :return:
        """
        def append_fifo(array, x):
            if len(array)<plotwin:
                array.append(x)
            else:
                array=np.hstack((array[1:], x))
            return array
        lines=[self.gluc, self.gluc_smooth, self.temp, self.times, self.refconc]
        lines=list(map(lambda x,y: append_fifo(x,y), lines, [gluc, gluc_smooth, temp, time, theoretical_conc]))

        self.templine.set_data(self.times, self.temp)
        self.glucline.set_data(self.times, self.gluc)
        self.theoretical_line.set_data(self.times, self.refconc)
        self.glucline_smooth.set_data(self.times, self.gluc_smooth)
        self.tempaxis.set_ylim((min(self.temp)-2, max(self.temp)+2))
        self.axis.set_ylim((min(self.gluc) - 2, max(self.gluc) + 2))
        if len(self.times)>5:
            self.axis.set_xlim((min(self.times), max(self.times)))
        else:
            self.axis.set_xlim((min(self.times)-timedelta(minutes=1), max(self.times)+timedelta(minutes=1)))
        #plt.draw()

    def terminate(self):
        plt.close('all')

    def call_back(self):
        '''called at every update '''
        while self.pipe.poll():
            command = self.pipe.recv()
            #print('received command',command)
            if command is None:
                self.terminate()
                return False
            else:
                self.update_plot(command['time'], command['gluc'], command['temp'], command['gluc_smooth'], command['theoretical_conc'])
        self.f.canvas.draw()
        return True

    def __call__(self, pipe):
        ''' starting point'''
        print('starting plotter...')
        self.pipe = pipe
        timer = self.f.canvas.new_timer(interval=100)
        timer.add_callback(self.call_back)
        timer.start()
        print('...done')
        plt.show()
        return


class MultiprocConnector(object):
    """ Connector between the plotter and the measurements. It allows real time updates"""
    def __init__(self):
        print ('multiproc created')
        self.plot_pipe, plotter_pipe = mp.Pipe()
        self.plotter = Plotter()
        self.plot_process = mp.Process(
            target=self.plotter, args=(plotter_pipe,), daemon=True)
        self.plot_process.start()

    def send_finished(self):
        self.plot_pipe.send(None)

    def update_data(self, time, gluc, temp, gluc_smooth, theoretical_conc):
        data = {'gluc': gluc, 'temp': temp, 'time': time, 'gluc_smooth': gluc_smooth, 'theoretical_conc':theoretical_conc}
        self.plot_pipe.send(data)
        return

class Measurement:
    """ Measurement object
        - initializes the serial communication with the glucose sensor
        - initializes the output csv file
        - handles addition of measurement points"""
    def __init__(self):
        self.calib=calibrationData()
        self.prevmeasurements=[]
        self.prev_measurement=None
        self.ser = self.init_serial()
        if saveOutput:
            self.outputfile = self.init_output()

    # def start_measuring(self):
    #     self.s = sched.scheduler(time.time, time.sleep)
    #     self.s.enter(Measurement_Rate, 1, self.append_measurement)
    #     self.s.run()

    def append_measurement(self):
        """
        Appending a measurement point:
            - resetting the serial input buffer (deleting the current one)
            - reading and processing 25 bytes
            - converting the data to glucose
            - applying a running median filter

        :return: timestamp, temperature, glucose with temperature compensation , smoothed glucose , theoretical concentration (if available)
        """
        ##empties the input buffer (no need to read the ones that are waiting)
        self.ser.reset_input_buffer()
        ## reads the next samples
        data = self.ser.read(25)
        timestamp = datetime.now()
        m = MeasurementPoint(data, timestamp)
        if m.validData:
            #calculating glucose
            m.gluc, m.gluc_comp = self.calib.calc_glucose(m)
            if len(self.prevmeasurements)<medianfiltwin:
                self.prevmeasurements.append(m.gluc_comp)
            else:
                self.prevmeasurements=np.hstack((self.prevmeasurements[1:], m.gluc_comp))
            m.gluc_smooth = np.median(self.prevmeasurements)
            #if self.prev_measurement is not None:
            #    m.gluc_smooth=self.prev_measurement*(1-lpf)+m.gluc_comp*lpf
            #else:
            #    m.gluc_smooth=m.gluc_comp
            #self.prev_measurement=m.gluc_smooth


            #self.measurements.append(m)
            print( timestamp, ' ch values:', m.channels,
                  'temperature:', m.temperature, 'glucose: {:.2f}'.format(m.gluc))


            m.theoretical_conc=np.nan
            if saveOutput and self.outputfile is not None:
                self.append_row(m)

        return timestamp, m.temperature, m.gluc_comp, m.gluc_smooth, m.theoretical_conc

    def init_output(self):
        """
        initializing the output csv file
        :return: the name of the output file
        """
        try:
            outputfile = datetime.strftime(datetime.now(), filePath + 'glucosesensor_%Y_%m_%d_%H_%M.csv')
            with open(outputfile, mode='w', newline='') as results:
                writer = csv.DictWriter(results, fieldnames=outputfields)
                writer.writeheader()
            return outputfile
        except:
            traceback.print_exc()
            return None

    def append_row(self, m):
        """
        Opens the output csv, append a row to it, and closes it
        :param m: MeasurementPoint
        :return: boolean whether it was successful
        """
        try:
            with open(self.outputfile, mode='a', newline='') as results:
                writer = csv.DictWriter(results, fieldnames=outputfields)
                towrite={'time': m.timestring,
                         'glucose [mM]': m.gluc_comp,
                         'glucose [mM] before temperature compensation':m.gluc,
                         'temperature [C]': m.temperature,
                         'temperature constant': self.calib.tempconstant
                         }
                writer.writerow(towrite)
                return True
        except:
            traceback.print_exc()
            return False

    def close(self):
        self.ser.close()
        print('closing serial connection')
        return

    def init_serial(self):
        """
        Initializing serial connection
        :return: serial object
        """
        try:
            if not fixedPort:
                portlist = [p.device for p in list(serial.tools.list_ports.comports()) if
                            'Prolific USB-to-Serial Comm Port' in p.description]
                port=portlist[0] if len(portlist)>0 else None
            else:
                port = comPort
            ser = serial.Serial(
                port=port,
                baudrate=9600,
                parity=serial.PARITY_NONE,
                stopbits=serial.STOPBITS_ONE,
                bytesize=serial.EIGHTBITS,
                timeout=1000)
            if ser.isOpen():
                print('closing previous session')
                ser.close()
                print('successfully closed')
            ser.open()
            print('successfully opened')
            return ser
        except:
            traceback.print_exc()
            exit()

def main():
    """main entry point"""
    try:
        connector = MultiprocConnector()
        #p=Plotter()
        mes=Measurement()
        if mes.ser is not None:
            print('starting measurements every {:.0f} seconds'.format(Measurement_Rate))
            while mes.ser.isOpen():
                t, temp, gluc, gluc_smooth, theoretical_conc = mes.append_measurement()
                connector.update_data(t, gluc, temp, gluc_smooth, theoretical_conc)
                time.sleep(Measurement_Rate)


    except KeyboardInterrupt:
        print('Interrupted')
        connector.send_finished()
        mes.close()
        exit()

    except:
        traceback.print_exc()
        exit()

if __name__ == '__main__':
    main()
