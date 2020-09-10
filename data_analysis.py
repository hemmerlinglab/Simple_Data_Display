#!/usr/local/bin/python3

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

		menu_bar = self.menuBar()
		menu_bar.setNativeMenuBar(False)
		file_menu = menu_bar.addMenu('File')
		file_menu.addAction(exit_action)
		controls_menu = menu_bar.addMenu('Controls')
		controls_menu.addAction(QAction('Up/Down = Move Time Cursors',self))
		controls_menu.addAction(QAction('Left/Right = Move Frequency Cursors',self))
		controls_menu.addAction(QAction('W/S = Increase/Decrease Time Cursor Difference',self))

		self.show()

	def closeEvent(self,event):
		reply = QuitMessage().exec_()
		if reply ==QMessageBox.Yes:
			event.accept()
		else:
			event.ignore()


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
		
		self.getButton = QPushButton('GET DATA')
		self.getButton.setFont(QFont('arial',fontsize))
		
		self.currTimeLab = QLabel('Curr Time:')
		self.currTimeLab.setFont(QFont('arial',fontsize))
		
		self.currTime = QLabel()
		self.currTime.setFont(QFont('arial',fontsize))
		
		self.currFreqLab = QLabel('Curr Freq:')
		self.currFreqLab.setFont(QFont('arial',fontsize))
		
		self.currFreq = QLabel()
		self.currFreq.setFont(QFont('arial',fontsize))

		GridLayout = QGridLayout()
		GridLayout.addWidget(self.baselab,0,0)
		GridLayout.addWidget(self.baseedit,0,1)
		GridLayout.addWidget(self.timelab,1,0)
		GridLayout.addWidget(self.timeedit,1,1)
		GridLayout.addWidget(self.getButton,2,0,1,2)
		GridLayout.addWidget(self.currTimeLab,3,0)
		GridLayout.addWidget(self.currTime,3,1)
		GridLayout.addWidget(self.currFreqLab,4,0)
		GridLayout.addWidget(self.currFreq,4,1)
		self.setLayout(GridLayout)


class CentralWidget(QWidget):
	def __init__(self,parent):
		QWidget.__init__(self,parent)
		self.basefilename = '20200616'
		self.time_stamp = '142554'
		self.filepath = self.basefilename
		self.f_idx = 20
		self.t_idx = 100
		self.t_off = 10
		self.t_avg_idx = self.t_idx + self.t_off
		self.offset = self.get_offset()
		self.setup()

	def setup(self):
		self.img_plot = MplCanvas(self,width=5,height=4,dpi=100)
		self.spec_plot = MplCanvas(self,width=5,height=1,dpi=100)
		self.shot_plot = MplCanvas(self,width=1,height=4,dpi=100)
		self.cont = ControlWidget(self)
		self.cont.getButton.clicked.connect(self.buttonClicked)
		self.cont.baseedit.setText(self.basefilename)
		self.cont.timeedit.setText(self.time_stamp)

		self.img_plot.keyPressEvent = self.keyPressEvent
		
		GridLayout = QGridLayout()
		GridLayout.addWidget(self.img_plot,0,0,3,6)
		GridLayout.addWidget(self.spec_plot,4,0,2,6)
		GridLayout.addWidget(self.shot_plot,0,6,3,2)
		GridLayout.addWidget(self.cont,4,6,2,2)
		self.setLayout(GridLayout)
		self.img_plot.setFocus()

	def make_img_plot(self):
		self.get_data()
		# plot_min = np.abs(self.ch0[0,0])
		self.img_plot.axes.cla()

		self.shot_min = np.min(self.ch0)
		self.shot_max = np.max(self.ch0[:,10])

		X,Y = np.meshgrid(self.freqs,self.times)
		self.img_plot.axes.pcolor(X,Y,self.ch0.T,vmin=self.shot_min,vmax=self.shot_max)
		self.img_plot.axes.set_title('{}_{}'.format(self.basefilename,self.time_stamp))
		
		self.vline = self.img_plot.axes.axvline(self.freqs[self.f_idx],alpha=0.4,color='red')
		self.vline2 = self.spec_plot.axes.axvline(self.freqs[self.f_idx],alpha=0.4,color='red')
		self.hline = self.img_plot.axes.axhline(self.times[self.t_idx],alpha=0.4,color='red')
		self.hline2 = self.img_plot.axes.axhline(self.times[self.t_avg_idx],alpha=0.4,color='red')
		self.hline3 = self.shot_plot.axes.axhline(self.times[self.t_idx],alpha=0.4,color='red')
		self.hline4 = self.shot_plot.axes.axhline(self.times[self.t_avg_idx],alpha=0.4,color='red')
		
		self.img_plot.axes.set_xlabel('Frequency (MHz from {:.6f})'.format(self.offset))
		self.img_plot.axes.set_ylabel('Time (ms)')
		
		self.img_plot.fig.canvas.draw()

	def make_slc_plots(self):
		self.t_avg_idx = self.t_idx + self.t_off
		self.spec_plot.axes.cla()
		self.shot_plot.axes.cla()

		self.spec_plot.axes.set_title('Spectrum')
		self.shot_plot.axes.set_xlabel('Signal (arb.)')
		# self.spec_plot.axes.set_ylabel('Time (ms)')
		self.shot_plot.axes.set_title('Single Shot')
		# self.shot_plot.axes.set_xlabel('Frequency (MHz from {:.6f})'.format(self.offset))
		self.spec_plot.axes.set_ylabel('Signal (arb.)')

		self.offset = self.get_offset()
		act_freq = self.offset*3 + (self.freqs[self.f_idx])*1e-6

		self.cont.currFreq.setText('{:.6f} THz'.format(act_freq))
		self.cont.currTime.setText('{} - {} ms'.format(str(self.times[self.t_idx]),str(self.times[self.t_avg_idx])))

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
		self.shot_plot.fig.canvas.draw()
		self.img_plot.fig.canvas.draw()

		self.repaint()

	def get_data(self):
		data_path = self.filepath
		time_stamp = self.time_stamp
		with open('{}/{}_{}_ch0_arr'.format(data_path,data_path,time_stamp),'r') as filename:
			ch0_data = np.genfromtxt(filename,delimiter=',')
		with open('{}/{}_{}_freqs'.format(data_path,data_path,time_stamp),'r') as filename:
			freqs = np.genfromtxt(filename)
		with open('{}/{}_{}_times'.format(data_path,data_path,time_stamp),'r') as filename:
			times = np.genfromtxt(filename)
		with open('{}/{}_{}_ch2_arr'.format(data_path,data_path,time_stamp),'r') as filename:
			ch2_data = np.genfromtxt(filename,delimiter=',')
		
		self.freqs = np.array(freqs)
		self.times = np.array(times)
		self.ch0 = self.get_avg(np.array(ch0_data),len(freqs))
		self.ch2 = self.get_avg(np.array(ch2_data),len(freqs))

		scl_fact = np.mean(self.ch0[:20,:])/np.mean(self.ch2[:20,:])

		self.ch0 -= scl_fact*self.ch2

		# for i in range(len(self.ch0[0,:])):
		# 	zero_offset = self.ch0[0,i]
		# 	self.ch0[:,i] -= zero_offset

	def get_offset(self):
		filename = '{}/{}_{}_conf'.format(self.filepath,self.basefilename,self.time_stamp)
		config = ConfigParser()
		config.read(filename)
		offset = np.float(config.get('offset_laser1','val')) #Thz
		return offset

	def get_avg(self,data,desired_len):
		new_data = np.zeros((desired_len,data.shape[1]))
		no_avgs = int(data.shape[0]/desired_len)
		for i in range(desired_len):
			for j in range(no_avgs):
				new_data[i,:] += data[(no_avgs*i+j),:]
		return new_data

	def buttonClicked(self):
		new_basefilename = self.cont.baseedit.text()
		new_timestamp = self.cont.timeedit.text()
		self.basefilename = new_basefilename
		self.time_stamp = new_timestamp
		self.filepath = new_basefilename
		self.get_data()
		self.make_img_plot()
		self.make_slc_plots()
		self.img_plot.setFocus()
		self.update()

	def keyPressEvent(self,event):
		if event.key() == Qt.Key_Up:
			self.t_idx += 1
			self.make_slc_plots()

		elif event.key() == Qt.Key_Down:
			self.t_idx -= 1
			self.make_slc_plots()

		elif event.key() == Qt.Key_Left:
			self.f_idx -= 1
			self.make_slc_plots()

		elif event.key() == Qt.Key_Right:
			self.f_idx += 1
			self.make_slc_plots()

		elif event.key() == Qt.Key_W:
			self.t_off += 1
			self.make_slc_plots()

		elif event.key() == Qt.Key_S:
			if self.t_off <= 1:
				pass
			else:
				self.t_off -= 1
			self.make_slc_plots()

		else:
			pass

		self.update()




if __name__ == '__main__':
	app = QApplication(sys.argv)
	main_window = MainWindow()
	app.exec_()