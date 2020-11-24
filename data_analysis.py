#!/usr/bin/python3

import sys
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
import matplotlib
matplotlib.use('QT5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import numpy as np
from configparser import ConfigParser
import time
from datetime import datetime

class MainWindow(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)
        self.setup()

    def setup(self):
        self.setGeometry(0,0,1400,800)
        self.setWindowTitle('Main Window')
        self.central_widget = CentralWidget(self)
        self.setCentralWidget(self.central_widget)

        exit_action = QAction('Quit',self)
        exit_action.triggered.connect(qApp.quit)

        grab_action = QAction('Save JPG',self)
        grab_action.triggered.connect(self.grab_screen)

        menu_bar = self.menuBar()
        menu_bar.setNativeMenuBar(False)
        file_menu = menu_bar.addMenu('File')
        file_menu.addAction(grab_action)
        file_menu.addAction(exit_action)
        
        controls_menu = menu_bar.addMenu('Spec Controls')
        controls_menu.addAction(QAction('Up/Down = Move Time Cursors',self))
        controls_menu.addAction(QAction('Left/Right = Move Frequency Cursors',self))
        controls_menu.addAction(QAction('W/S = Increase/Decrease Time Cursor Difference',self))

        target_menu = menu_bar.addMenu('Target Controls')
        target_menu.addAction(QAction('Arrows = Move X/Y Position',self))
        target_menu.addAction(QAction('W/S = Increase/Decrease Time Cursor Difference',self))
        target_menu.addAction(QAction('A/D = Move Time Cursors (invisible)',self))

        self.show()

    def closeEvent(self,event):
        reply = QuitMessage().exec_()
        if reply ==QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()

    def grab_screen(self):
        date = datetime.now()
        nowtime = date.strftime('%Y%m%d_%H%M%S:')
        filename = 'Screenshots/{}{}_{}_{}_{}_{}.jpg'.format(nowtime,self.central_widget.basefilename,self.central_widget.time_stamp,str(self.central_widget.f_idx),str(self.central_widget.t_idx),str(self.central_widget.t_off))
        screen = QApplication.primaryScreen()
        screenshot = screen.grabWindow(self.central_widget.winId())
        screenshot.save(filename,'jpg')


class QuitMessage(QMessageBox):
    def __init__(self):
        QMessageBox.__init__(self)
        self.setText('Are you sure you want to quit?')
        self.addButton(self.No)
        self.addButton(self.Yes)


class MplCanvas(FigureCanvasQTAgg):
    def __init__(self,parent=None,width=5,height=4,dpi=100):
        self.fig = Figure(figsize=(width,height),dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        super(MplCanvas,self).__init__(self.fig)


class BlinkButton(QPushButton):
    def __init__(self, *args, **kwargs):
        QPushButton.__init__(self, *args, **kwargs)
        self.default_color = self.getColor()

    def getColor(self):
        return self.palette().color(QPalette.Button)

    def setColor(self, value):
        if value == self.getColor():
            return
        palette = self.palette()
        palette.setColor(self.backgroundRole(), value)
        self.setAutoFillBackground(True)
        self.setPalette(palette)

    def reset_color(self):
        self.setColor(self.default_color)

    color = pyqtProperty(QColor, getColor, setColor)


class ControlWidget(QWidget):
    def __init__(self,parent):
        QWidget.__init__(self,parent)
        self.setup()

    def setup(self):
        fontsize = 20
        self.baselab = QLabel('Base File Name:')
        self.baselab.setFont(QFont('arial',fontsize))

        self.baseedit = QLineEdit()
        self.baseedit.setFont(QFont('arial',fontsize))
        
        self.timelab = QLabel('Timestamp:')
        self.timelab.setFont(QFont('arial',fontsize))
        
        self.timeedit = QLineEdit()
        self.timeedit.setFont(QFont('arial',fontsize))
        
        self.getButton = BlinkButton('GET DATA')
        self.getButton.setFont(QFont('arial',fontsize))
        
        self.currTimeLab = QLabel('Current Time:')
        self.currTimeLab.setFont(QFont('arial',fontsize))
        
        self.currTime = QLabel()
        self.currTime.setFont(QFont('arial',fontsize))
        
        self.currFreqLab = QLabel('Current Freq:')
        self.currFreqLab.setFont(QFont('arial',fontsize))
        
        self.currFreq = QLabel()
        self.currFreq.setFont(QFont('arial',fontsize))

        self.multFactLab = QLabel('Frequency Multiple:')
        self.multFactLab.setFont(QFont('arial',fontsize))

        self.multFact = QLineEdit()
        self.multFact.setFont(QFont('arial',fontsize))

        GridLayout = QGridLayout()
        GridLayout.addWidget(self.baselab,0,0)
        GridLayout.addWidget(self.baseedit,0,1)
        GridLayout.addWidget(self.timelab,1,0)
        GridLayout.addWidget(self.timeedit,1,1)
        GridLayout.addWidget(self.multFactLab,2,0)
        GridLayout.addWidget(self.multFact,2,1)
        GridLayout.addWidget(self.getButton,3,0,1,2)
        GridLayout.addWidget(self.currTimeLab,4,0)
        GridLayout.addWidget(self.currTime,4,1)
        GridLayout.addWidget(self.currFreqLab,5,0)
        GridLayout.addWidget(self.currFreq,5,1)
        self.setLayout(GridLayout)


class CentralWidget(QWidget):
    def __init__(self,parent):
        QWidget.__init__(self,parent)
        self.basefilename = '20201030'
        self.time_stamp = '151550'
        self.freq_mult = 1
        self.filepath = '/home/molecules/software/data/'+self.basefilename
        # self.filepath = self.basefilename
        self.f_idx = 10
        self.t_idx = 100
        self.t_off = 20
        self.t_avg_idx = self.t_idx + self.t_off
        self.set_filetype()
        self.offset = self.get_offset()
        self.setup()
        self.default_button_color = self.palette().color(QPalette.Button)

    def setup(self):
        self.img_plot = MplCanvas(self,width=5,height=4,dpi=100)
        self.spec_plot = MplCanvas(self,width=5,height=1,dpi=100)
        self.shot_plot = MplCanvas(self,width=1,height=4,dpi=100)
        self.cont = ControlWidget(self)
        self.cont.getButton.clicked.connect(self.buttonClicked)
        self.cont.baseedit.setText(self.basefilename)
        self.cont.timeedit.setText(self.time_stamp)
        self.cont.multFact.setText(str(self.freq_mult))

        self.img_plot.keyPressEvent = self.keyPressEvent
        
        GridLayout = QGridLayout()
        GridLayout.addWidget(self.img_plot,0,0,3,6)
        GridLayout.addWidget(self.spec_plot,4,0,2,6)
        GridLayout.addWidget(self.shot_plot,0,6,3,2)
        GridLayout.addWidget(self.cont,4,6,2,2)
        self.setLayout(GridLayout)
        self.img_plot.setFocus()
        self.update()

    def make_img_plot(self):
        # self.get_data()
        # plot_min = np.abs(self.ch0[0,0])
        self.img_plot.axes.cla()

        self.shot_min = np.min(self.ch0)
        self.shot_max = np.max(self.ch0[:,10])

        if self.filetype == 'spec':
            X,Y = np.meshgrid(self.freqs,self.times)
            self.img_plot.axes.pcolor(X,Y,self.ch0.T,vmin=self.shot_min,vmax=self.shot_max)
            self.img_plot.axes.set_title('{}_{} Spectrum'.format(self.basefilename,self.time_stamp))

            self.vline = self.img_plot.axes.axvline(self.freqs[self.f_idx],alpha=0.4,color='red')
            self.vline2 = self.spec_plot.axes.axvline(self.freqs[self.f_idx],alpha=0.4,color='red')
            self.hline = self.img_plot.axes.axhline(self.times[self.t_idx],alpha=0.4,color='red')
            self.hline2 = self.img_plot.axes.axhline(self.times[self.t_avg_idx],alpha=0.4,color='red')
            self.hline3 = self.shot_plot.axes.axhline(self.times[self.t_idx],alpha=0.4,color='red')
            self.hline4 = self.shot_plot.axes.axhline(self.times[self.t_avg_idx],alpha=0.4,color='red')
            
            self.img_plot.axes.set_xlabel('Frequency (MHz from {:.6f})'.format(self.offset*self.freq_mult))
            self.img_plot.axes.set_ylabel('Time (ms)')

        elif self.filetype == 'target':
            tar_img = self.get_tar_img()
            self.img_plot.axes.pcolor(self.xpos,self.ypos,tar_img.T,vmin=self.shot_min,vmax=self.shot_max)
            self.img_plot.axes.set_title('{}_{} Target'.format(self.basefilename,self.time_stamp))

            self.vline = self.img_plot.axes.axvline(self.xpos[0,self.x_idx],alpha=0.4,color='red')
            self.vline2 = self.img_plot.axes.axvline(self.xpos[0,self.x_idx+1],alpha=0.4,color='red')
            self.hline = self.img_plot.axes.axhline(self.ypos[self.y_idx,0],alpha=0.4,color='red')
            self.hline2 =  self.img_plot.axes.axhline(self.ypos[self.y_idx+1,0],alpha=0.4,color='red')
            self.hline3 = self.shot_plot.axes.axhline(self.times[self.t_idx],alpha=0.4,color='red')
            self.hline4 = self.shot_plot.axes.axhline(self.times[self.t_avg_idx],alpha=0.4,color='red')

            self.img_plot.axes.invert_yaxis()

            self.img_plot.axes.set_xlabel('X Position')
            self.img_plot.axes.set_ylabel('Y Position')

        else:
            print('AAAAAHHHH')
        
        
        self.img_plot.fig.canvas.draw()
        self.repaint()

    def get_tar_img(self):
        # print(self.ch0.shape)
        slc_arr = np.mean(self.ch0[:,self.t_idx:self.t_avg_idx],axis=1)
        # print(slc_arr.shape)
        return slc_arr.reshape((self.steps_x,self.steps_y))

    def get_tar_slc(self):
        new_arr = np.zeros((self.steps_x,self.steps_y,self.ch0.shape[1]))
        for i in range(self.steps_x):
            for j in range(self.steps_y):
                new_arr[i,j,:] = self.ch0[j+i*self.steps_y,:]
        return new_arr

    def make_slc_plots(self):
        self.t_avg_idx = self.t_idx + self.t_off
        self.shot_plot.axes.cla()

        if self.filetype == 'spec':
            self.spec_plot.axes.cla()
            self.spec_plot.axes.set_title('Spectrum')
            self.shot_plot.axes.set_xlabel('Signal (arb.)')
            self.shot_plot.axes.set_title('Single Shot')

        # self.spec_plot.axes.set_ylabel('Time (ms)')
        # self.shot_plot.axes.set_xlabel('Frequency (MHz from {:.6f})'.format(self.offset))
            self.spec_plot.axes.set_ylabel('Signal (arb.)')

            act_freq = (self.offset + (self.freqs[self.f_idx])*1e-6)*self.freq_mult

            self.cont.currFreq.setText('{:.6f} THz'.format(act_freq))
            self.cont.currTime.setText('{:.3f} - {:.3f} ms'.format(self.times[self.t_idx],self.times[self.t_avg_idx]))

            self.spec_plot.axes.plot(self.freqs,np.mean(self.ch0[:,self.t_idx:self.t_avg_idx],axis=1))
            self.shot_plot.axes.plot(self.ch0[self.f_idx,:],self.times)
            self.spec_plot.axes.set_xlim(np.min(self.freqs),np.max(self.freqs))
            self.shot_plot.axes.set_xlim(self.shot_min,self.shot_max)
            self.shot_plot.axes.set_ylim(np.min(self.times),np.max(self.times))
            self.spec_plot.axes.set_ylim(self.shot_min,self.shot_max)
            self.vline.remove()
            self.vline2.remove()
            self.hline.remove()
            self.hline2.remove()
            self.hline3.remove()
            self.hline4.remove()

            self.vline = self.img_plot.axes.axvline(self.freqs[self.f_idx],alpha=0.4,color='red')
            self.vline2 = self.spec_plot.axes.axvline(self.freqs[self.f_idx],alpha=0.4,color='red')
            self.hline = self.img_plot.axes.axhline(self.times[self.t_idx],alpha=0.4,color='red')
            self.hline2 = self.img_plot.axes.axhline(self.times[self.t_avg_idx],alpha=0.4,color='red')
            self.hline3 = self.shot_plot.axes.axhline(self.times[self.t_idx],alpha=0.4,color='red')
            self.hline4 = self.shot_plot.axes.axhline(self.times[self.t_avg_idx],alpha=0.4,color='red')
        
            self.spec_plot.fig.canvas.draw()

        elif self.filetype == 'target':
            self.shot_plot.axes.set_title('X={:.2f}/Y={:.2f} Shot'.format(self.xpos[0,self.x_idx],self.ypos[self.y_idx,0]))

            slc_plot = self.get_tar_slc()
            self.shot_plot.axes.plot(slc_plot[self.x_idx,self.y_idx,:],self.times)
            self.shot_plot.axes.set_xlim(self.shot_min,self.shot_max)
            self.shot_plot.axes.set_ylim(np.min(self.times),np.max(self.times))

            self.vline.remove()
            self.vline2.remove()
            self.hline.remove()
            self.hline2.remove()
            self.hline3.remove()
            self.hline4.remove()
            self.vline = self.img_plot.axes.axvline(x=self.xpos[0,self.x_idx],alpha=0.4,color='red')
            self.vline2 = self.img_plot.axes.axvline(x=self.xpos[0,self.x_idx+1],alpha=0.4,color='red')
            self.hline = self.img_plot.axes.axhline(y=self.ypos[self.y_idx,0],alpha=0.4,color='red')
            self.hline2 =  self.img_plot.axes.axhline(y=self.ypos[self.y_idx+1,0],alpha=0.4,color='red')
            self.hline3 = self.shot_plot.axes.axhline(y=self.times[self.t_idx],alpha=0.4,color='red')
            self.hline4 = self.shot_plot.axes.axhline(y=self.times[self.t_avg_idx],alpha=0.4,color='red')

        else:
            print('OOOOF!')

        self.shot_plot.fig.canvas.draw()
        self.img_plot.fig.canvas.draw()

        self.repaint()

    def get_data(self):
        self.offset = self.get_offset()
        #data_path = self.basefilename
        #time_stamp = self.time_stamp
        with open('{}/{}_{}_ch0_arr'.format(self.filepath,self.basefilename,self.time_stamp),'r') as filename:
            ch0_data = np.genfromtxt(filename,delimiter=',')
        with open('{}/{}_{}_times'.format(self.filepath,self.basefilename,self.time_stamp),'r') as filename:
            times = np.genfromtxt(filename)
        with open('{}/{}_{}_ch3_arr'.format(self.filepath,self.basefilename,self.time_stamp),'r') as filename:
            ch3_data = np.genfromtxt(filename,delimiter=',')

        if self.filetype == 'target':
            with open('{}/{}_{}_posx'.format(self.filepath,self.basefilename,self.time_stamp),'r') as filename:
                x_data = np.genfromtxt(filename,delimiter=',')
            with open('{}/{}_{}_posy'.format(self.filepath,self.basefilename,self.time_stamp),'r') as filename:
                y_data = np.genfromtxt(filename,delimiter=',')

            self.xpos = np.array(x_data)
            self.ypos = np.array(y_data)
            self.xpos = self.xpos.reshape((self.steps_x,self.steps_y))
            self.ypos = self.ypos.reshape((self.steps_x,self.steps_y))
            self.x_idx = int(self.steps_x/2)
            self.y_idx = int(self.steps_y/2)
            # print(self.xpos.shape)
            # print(self.ypos.shape)

        if self.filetype == 'spec':
            with open('{}/{}_{}_freqs'.format(self.filepath,self.basefilename,self.time_stamp),'r') as filename:
                freqs = np.genfromtxt(filename)
            self.freqs = np.array(freqs)*self.freq_mult
            # print(self.freqs.shape)
        
        self.times = np.array(times)

        # print(ch0_data.shape)

        self.ch0 = self.get_avg(np.array(ch0_data),int(len(ch0_data)/self.no_avgs))
        self.ch3 = self.get_avg(np.array(ch3_data),int(len(ch0_data)/self.no_avgs))

        # print(self.ch0.shape)

        scl_fact = np.mean(self.ch0[:20,:])/np.mean(self.ch3[:20,:])

        #self.ch0 -= scl_fact*self.ch3

        self.ch0 = self.rem_offset(self.ch0)

        self.f_idx = int(self.ch0.shape[0]/2)
    
        # print(self.filetype)

        # print(self.ch0.shape)
        # print(self.times.shape)

        # for i in range(len(self.ch0[0,:])):
        #   zero_offset = self.ch0[0,i]
        #   self.ch0[:,i] -= zero_offset

    def set_filetype(self):
        filename = '{}/{}_{}_conf'.format(self.filepath,self.basefilename,self.time_stamp)
        config = ConfigParser()
        config.read(filename)
        try:
            self.min_x = np.float(config.get('min_x','val'))
            self.filetype = 'target'
        except:
            self.filetype = 'spec'

    def get_offset(self):
        filename = '{}/{}_{}_conf'.format(self.filepath,self.basefilename,self.time_stamp)
        config = ConfigParser()
        config.read(filename)
        self.no_avgs = int(config.get('scan_count','val'))
        if self.filetype == 'spec':
            self.whichlaser = int(config.get('which_scanning_laser','val'))
            offset = np.float(config.get('offset_laser1','val'))
            if self.whichlaser == 2:
                offset = np.float(config.get('offset_laser2','val'))
        elif self.filetype == 'target':
            offset = np.float(config.get('cooling_set','val'))
            self.steps_x = int(float(config.get('steps_x','val')))
            self.steps_y = int(float(config.get('steps_y','val')))
        else:
            print('AAAAAAAAAAAAAAAAAAAAAAAAAA')
        return offset

    def get_avg(self,data,desired_len):
        new_data = np.zeros((desired_len,data.shape[1]))
        no_avgs = int(data.shape[0]/desired_len)
        for i in range(desired_len):
            for j in range(no_avgs):
                new_data[i,:] += self.rem_offset(data[(no_avgs*i+j),:])
        return new_data

    def rem_offset(self,data):
        new_data = np.zeros(data.shape)
        try:
            for i in range(new_data.shape[0]):
                new_data[i,:] = data[i,:] - np.mean(data[i,0:20])
        except:
            new_data = data - np.mean(np.array([np.mean(data[0:20]),np.mean(data[0:-20])]))
        return new_data

    def buttonClicked(self):
        try:
            self.cont.getButton.reset_color()
            self.cont.getButton.setText('GET DATA')
            self.cont.getButton.setFont(QFont('arial',20))
            new_basefilename = self.cont.baseedit.text()
            new_timestamp = self.cont.timeedit.text()
            new_freq_mult = int(self.cont.multFact.text())
            self.basefilename = new_basefilename
            self.time_stamp = new_timestamp
            self.freq_mult = new_freq_mult
            self.filepath = '/home/molecules/software/data/'+new_basefilename
            # self.filepath = new_basefilename
            # self.get_data()
            self.set_filetype()
            self.get_data()
            self.update()
            self.make_img_plot()
            self.make_slc_plots()
            self.img_plot.setFocus()
        except:
            self.cont.getButton.setColor(Qt.red)
            self.cont.getButton.setText('ERROR')
            self.cont.getButton.repaint()
            self.cont.getButton.update()

        self.update()


    def keyPressEvent(self,event):
        if event.key() == Qt.Key_Up:
            if self.filetype == 'spec':
                if self.t_idx+self.t_off < self.ch0.shape[1]-1:
                    self.t_idx += 1
            elif self.filetype == 'target':
                if self.y_idx > 0:
                    self.y_idx -= 1
            else:
                pass
            self.make_slc_plots()

        elif event.key() == Qt.Key_Down:
            if self.filetype == 'spec':
                if self.t_idx > 0:
                    self.t_idx -= 1
            elif self.filetype == 'target':
                if self.y_idx+1 < self.steps_y-1:
                    self.y_idx += 1
            else:
                pass
            self.make_slc_plots()

        elif event.key() == Qt.Key_Left:
            if self.filetype == 'spec':
                if self.f_idx > 0:
                    self.f_idx -= 1
            elif self.filetype == 'target':
                if self.x_idx > 0:
                    self.x_idx -= 1
            else:
                pass
            self.make_slc_plots()

        elif event.key() == Qt.Key_Right:
            if self.filetype == 'spec':
                if self.f_idx < self.ch0.shape[0]-1:
                    self.f_idx += 1
            elif self.filetype == 'target':
                if self.x_idx < self.steps_x-1:
                    self.x_idx += 1
            else:
                pass
            self.make_slc_plots()

        elif event.key() == Qt.Key_W:
            if self.t_idx+self.t_off < self.ch0.shape[1]:
                self.t_off += 1
                self.make_slc_plots()
                if self.filetype == 'target':
                    self.make_img_plot()

        elif event.key() == Qt.Key_S:
            if self.t_off <= 1:
                pass
            else:
                self.t_off -= 1
                self.make_slc_plots()
                if self.filetype == 'target':
                    self.make_img_plot()

        elif event.key() == Qt.Key_A:
            if self.filetype == 'target':
                if self.t_idx > 0:
                    self.t_idx -= 1
                self.make_slc_plots()
                self.make_img_plot()

        elif event.key() == Qt.Key_D:
            if self.filetype == 'target':
                if self.t_idx+self.t_off < self.ch0.shape[1]-1:
                    self.t_idx += 1
                self.make_slc_plots()
                self.make_img_plot()

        else:
            pass

        self.update()




if __name__ == '__main__':
    app = QApplication(sys.argv)
    main_window = MainWindow()
    app.exec_()
