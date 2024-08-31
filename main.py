# importing the required module  
import sys
import time
import numpy as np
import pandas as pd
from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt, QSize, QRect, QPoint, QPointF
from PyQt5.QtWidgets import QApplication, QLabel, QSizePolicy, QScrollArea, QMessageBox, QMainWindow, QMenu, QAction, qApp, QFileDialog, QToolBar
from PyQt5.QtGui import QImage, QMouseEvent, QPixmap, QPalette, QPainter, QIcon, QKeySequence, QPen
from PyQt5.QtPrintSupport import QPrintDialog, QPrinter
from scipy.spatial import KDTree

class QImageViewer(QMainWindow):
    # defining the initializing function
    def __init__(self):
        super().__init__()
        # configuring the window
        self.setWindowTitle('Pratyaksa - Spatial Data Analysis Tool')
        # configuring the width and height of the window  
        self.window_width, self.window_height = self.geometry().width(), self.geometry().height()
        # setting the icon of the window
        self.setWindowIcon(QIcon('./utils/icon.png'))
        # using the resize() to set the size of the application  
        self.resize(self.window_width * 2, self.window_height * 2)
        # creating an object of the QPrinter class
        self.printerObj = QPrinter()
        # setting the initial scale factor
        self.scale_factor = 0.0
        # creating an object of the QLabel class to display the label
        self.image_label = QLabel()
        # setting the background color of the label to display the image using the setBackgoundRole() method and QPalette class
        self.image_label.setBackgroundRole(QPalette.Base)
        # setting the size policy of the label using the setSizePolicy() method and QSizePolicy class
        self.image_label.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Ignored)
        #setting the setScaledContents() method to True to manually adjust the aspect ratio of the image in the application
        self.image_label.setScaledContents(True)
        # creating an object of the QScrollArea class to display the scroll area
        self.scroll_area = QScrollArea()
        # setting the background color of the scroll to discplay the background using the setBackgoundRole() method and QPalette class 
        self.scroll_area.setBackgroundRole(QPalette.Dark)
        # setting the scrolling area to the image label using the setWidget() method  
        self.scroll_area.setWidget(self.image_label) 
        # setting the visibility of the scrolling area with the help of the setVisible() method  
        self.scroll_area.setVisible(False)
        # setting the central widget to the scroll area using the setCentralWidget() method 
        self.setCentralWidget(self.scroll_area)
        # atrributes for dragging and hand tool functionality
        self.dragging = False
        self.lastMousePosition = None
        # Make the image label accept mouse events
        self.image_label.setMouseTracking(True)
        self.setMouseTracking(True)
        # to store the path of the loaded csv and the tissue possitions from the csv
        self.csvFilePath = None
        self.tissuePositions = None
        self.tissue_positions_complete=None
        self.isDrawing = False
        self.drawingToolActive = False
        self.originalPixmap = None  
        self.drawablePixmap = None
        self.highlightLayer = None # for the backtrack tool
        self.drawnPixels = pd.DataFrame(columns=['x', 'y'])
        # Initialize list to store drawn points
        self.drawnPoints = []

        #--------------------------------------  
        # Creating a File Menu  
        #--------------------------------------  
        self.filemenu = self.menuBar().addMenu('&File')  
        #--------------------------------------  
        # Creating a File Toolbar  
        #--------------------------------------  
        self.filetoolbar = QToolBar('File')  
        self.filetoolbar.setIconSize(QSize(30, 30))  
        self.addToolBar(Qt.ToolBarArea.TopToolBarArea, self.filetoolbar)
        # creating the menu options
        # calling the user defined makeAction() method to create the action for the menu options
        self.open_doc_opt = self.makeAction(self, 'utils/openImage.ico', 'Open Image...', 'Open Image...', self.openImage)
        # using the setShortcut() method to set a shortcut to execute the 'Open' command  
        self.open_doc_opt.setShortcut(QKeySequence.Open)
        # calling the user-defined makeAction() method to create the action to print the file  
        self.print_opt = self.makeAction(self, './utils/printer.ico', 'Print', 'Print', self.printImage)
        #using the setShortcut() method to set a shortcut to execute the 'Print' command
        self.print_opt.setShortcut(QKeySequence.Print)
        # initially disabling the action by setting the value of setEnabled() method to False  
        self.print_opt.setEnabled(False)
        # using the addActions() method to add all the created actions to the 'File' menu and toolbar 
        self.filemenu.addActions([self.open_doc_opt, self.print_opt])
        self.filetoolbar.addActions([self.open_doc_opt, self.print_opt])  
        # adding the separator  
        self.filemenu.addSeparator()  
        # calling the user-defined makeAction() method to create the action to close the application  
        self.exit_opt = self.makeAction(self, '', 'Exit', 'Exit', self.close) 
        # using the setShortcut() method to set a shortcut to execute the 'Close' command  
        self.print_opt.setShortcut(QKeySequence.Close)
        # using the addActions() method to add all the created actions to the 'File' menu and toolbar  
        self.filemenu.addActions([self.exit_opt])
        #setting up extract tool
        self.extractToolAction = self.makeAction(self, './utils/extract.ico', 'Extract', 'Extract data from image', self.extractBarcodes)
        self.filemenu.addAction(self.extractToolAction)
        self.filetoolbar.addAction(self.extractToolAction)
        #--------------------------------------  
        # Creating a View Menu  
        #--------------------------------------  
        self.viewmenu = self.menuBar().addMenu('&View')  
        #--------------------------------------  
        # Creating an View Tool bar  
        #--------------------------------------  
        self.viewtoolbar = QToolBar('Edit')  
        self.viewtoolbar.setIconSize(QSize(30, 30))  
        self.addToolBar(Qt.ToolBarArea.TopToolBarArea, self.viewtoolbar) 
        self.zoomIN_opt = self.makeAction(self, './utils/zoomIn.ico', 'Zoom In (25%)', 'Zoom In (25%)', self.zoom_in)
        # calling the user-defined makeAction() method to create the action to zoom in the image
        self.zoomIN_opt.setShortcut(QKeySequence.ZoomIn)
        # initially disabling the action by setting the value of setEnabled() method to False
        self.zoomIN_opt.setEnabled(False)
        # calling the user-defined makeAction() method to create the action to zoom out the image  
        self.zoomOUT_opt = self.makeAction(self, './utils/zoomOut.ico', 'Zoom Out (25%)', 'Zoom Out (25%)', self.zoom_out) 
        # using the setShortcut() method to set a shortcut to execute the 'Zoom Out' command  
        self.zoomOUT_opt.setShortcut(QKeySequence.ZoomOut)  
        # initially disabling the action by setting the value of setEnabled() method to False  
        self.zoomOUT_opt.setEnabled(False)  
        # calling the user-defined makeAction() method to create the action to set the normal size of the image  
        self.normalSize_opt = self.makeAction(self, '', 'Normal Size', 'Normal Size', self.normal_size)
        # using the setShortcut() method to set a shortcut to execute the 'Normal Size' command
        self.normalSize_opt.setShortcut(QKeySequence('Ctrl+S'))
        # initially disabling the action by setting the value of setEnabled() method to False
        self.normalSize_opt.setEnabled(False)
        # setting the initial value of setCheckable() method to True
        self.normalSize_opt.setCheckable(True)
        # using the addActions() method to add all the created actions to the 'View' menu and toolbar  
        self.viewmenu.addActions([self.zoomIN_opt, self.zoomOUT_opt, self.normalSize_opt])
        self.viewtoolbar.addActions([self.zoomIN_opt, self.zoomOUT_opt])  
        # adding the separator  
        self.viewmenu.addSeparator()  
        self.viewtoolbar.addSeparator()  
        # calling the user-defined makeAction() method to create the action to fit the image to the window  
        self.fitToWindow_opt = self.makeAction(self, './utils/fitToWindow.ico', 'Fit To Window', 'Fit To Window', self.fit_to_window)
        # using the setShortcut() method to set a shortcut to execute the 'Fit To Window' command
        self.fitToWindow_opt.setShortcut(QKeySequence('Ctrl+F'))
        # initially disabling the action by setting the value of setEnabled() method to False
        self.fitToWindow_opt.setEnabled(False)
        # using the addActions() method to add all the created actions to the 'View' menu and toolbar
        self.viewmenu.addActions([self.fitToWindow_opt])
        self.viewtoolbar.addActions([self.fitToWindow_opt])
        self.lastUpdateTime = None

        # Dictionary to store loaded images for caching
        self.image_cache = {}

        #--------------------------------------  
        # Creating a Tools Menu  
        #--------------------------------------  
        self.toolsmenu = self.menuBar().addMenu('&Tools')  
        #--------------------------------------  
        # Creating an View Tool bar  
        #--------------------------------------      
        self.toolstoolbar = QToolBar('Tools')
        self.toolstoolbar.setIconSize(QSize(30, 30))
        self.addToolBar(Qt.ToolBarArea.TopToolBarArea, self.toolstoolbar)
        # setting up tools
        self.setupTools()
        
    # defining the required methods of the class
    # defining the method to open the image file
    def openImage(self):
        # creating an object of the QFileDialog.Options class  
        selections = QFileDialog.Options()  
        # calling the getOpenFileName() method to browse the image from the directory  
        file_name, _ = QFileDialog.getOpenFileName(  
            self,  
            'QFileDialog.getOpenFileName()',  
            '',  
            'Images (*.png *.jpeg *.jpg *.bmp *.gif)',  
            options = selections  
            )  
        # if the file name is not an empty string  
        if file_name:  
            pixmap = self.loadImage(file_name)
            if pixmap:
                self.image_label.setPixmap(pixmap)
                # enabling the visibility of the scroll area  
                self.scroll_area.setVisible(True)  
                # enabling the "Print" action  
                self.print_opt.setEnabled(True)  
                # calling the fit_to_window() method  
                self.fit_to_window()  
                # enabling the "Fit To Window" action  
                self.fitToWindow_opt.setEnabled(True)  
                # calling the update_actions() method  
                self.update_actions() 
                # enabling the "Draw" action 
                self.drawToolAction.setEnabled(True)

                # if the "Fit To Window" action is not checked  
                if not self.fitToWindow_opt.isChecked():  
                    # calling the adjustSize() method to adjust the size of the image  
                    self.image_label.adjustSize() 

    def loadImage(self, file_name):
        if file_name in self.image_cache:
            self.originalPixmap = self.image_cache[file_name]
        else:
            image = QImage(file_name)
            if image.isNull():
                QMessageBox.information(self, 'Image Viewer', 'Cannot load %s.' % file_name)
                return None
            self.originalPixmap = QPixmap.fromImage(image)
            self.image_cache[file_name] = self.originalPixmap

        self.drawablePixmap = QPixmap(self.originalPixmap.size())
        self.drawablePixmap.fill(Qt.transparent)  # Use transparent for drawing on top

        self.updateImageDisplay()
        return self.drawablePixmap and self.originalPixmap

    def updateImageDisplay(self):
        # Ensure the original pixmap is not empty
        if not self.originalPixmap:
            return
        
        # Create a combined pixmap to merge the layers
        combinedPixmap = QPixmap(self.originalPixmap.size())
        combinedPixmap.fill(Qt.transparent)
        
        # Merge the original, drawable, and highlight layers
        painter = QPainter(combinedPixmap)
        painter.drawPixmap(0, 0, self.originalPixmap)
        painter.drawPixmap(0, 0, self.drawablePixmap)
        if self.highlightLayer:
            painter.drawPixmap(0, 0, self.highlightLayer)
        painter.end()
        
        self.image_label.setPixmap(combinedPixmap)


    #defining the method to print the image
    def printImage(self):
        # creating an object of the QPrintDialog class to print the image  
        print_dialog = QPrintDialog(self.printerObj, self)  
        # if the dialog is accepted  
        if print_dialog.exec_():  
            the_painter = QPainter(self.printerObj)  
            # creating a rectangle to place the image  
            rectangle = the_painter.viewport()  
            # defining the size of the image  
            the_size = self.image_label.pixmap().size()  
            # scaling the image to the Aspect Ratio  
            the_size.scale(rectangle.size(), Qt.KeepAspectRatio)  
            # setting the viewport of the image by calling the setViewport() method
            the_painter.setViewport(rectangle.x(), rectangle.y(), the_size.width(), the_size.height())  
            # calling the setWindow() method  
            the_painter.setWindow(self.image_label.pixmap().rect())  
            # calling the drawPixmap() method  
            the_painter.drawPixmap(0, 0, self.image_label.pixmap())

    #defining the zoom in on the image
    def zoom_in(self):
        # calling the scale_image() method to scale the image  
        self.scale_image(1.25)

    #defining the zoom out on the image
    def zoom_out(self):
        # calling the scale_image() method to scale the image  
        self.scale_image(0.8)

    #defining the method to set the normal size of the image
    def normal_size(self):
        #calling the adjustSize() method to adjust the size of the image
        self.image_label.adjustSize()
        #setting the scale factor to 1.0
        self.scale_factor = 1.0

    #defining the method to fit the image to the window
    def fit_to_window(self):
        #retriving the boolean value from the "Fit to window" action
        fitToWindow = self.fitToWindow_opt.isChecked()
        #configuring the scroll area to resizable
        self.scroll_area.setWidgetResizable(fitToWindow)
        #if the retrived value is False, calling the user-defined normal_size() method
        if not fitToWindow:
            self.normal_size()
        # calling the user defined update_actions() method
        self.update_actions()

    #defining the method to scale the image  
    def update_actions(self):
        # enabling the "Zoom In", "Zoom Out", and "Normal Size" actions, if the "Fir To Window" is unchecked
        self.zoomIN_opt.setEnabled(not self.fitToWindow_opt.isChecked()) 
        self.zoomOUT_opt.setEnabled(not self.fitToWindow_opt.isChecked())  
        self.normalSize_opt.setEnabled(not self.fitToWindow_opt.isChecked())   

    # defining the method to scale the image  
    # defining the method to scale the image  
    def scale_image(self, sf):
        # defining the scaling factor of the image  
        self.scale_factor *= sf  
        # using the resize() method to resize the image as per the scaling factor 
        self.image_label.resize(self.scale_factor * self.image_label.pixmap().size())  
        # calling the user-defined adjust_scroll_bar() method to adjust the scrollbar as per the scaling factor  
        self.adjust_scroll_bar(self.scroll_area.horizontalScrollBar(), sf)
        self.adjust_scroll_bar(self.scroll_area.verticalScrollBar(), sf)  
        # toggling the "Zoom In" and "Zoom Out" actions without strict limits
        self.zoomIN_opt.setEnabled(True)  # Always allow zooming in
        self.zoomOUT_opt.setEnabled(True)  # Remove lower limit for zooming out


    #defining the method to adjust the scrollbar
    def adjust_scroll_bar(self, scroll_bar, sf):  
        # setting the value of the scrollbar to the minimum value  
        scroll_bar.setValue(int(sf * scroll_bar.value() + ((sf - 1) * scroll_bar.pageStep() / 2)))
    
    def mousePressEvent(self, event):
        super().mousePressEvent(event)  # Call the parent class's method
        if event.button() == Qt.LeftButton and self.drawingToolActive:
            self.isDrawing = True
            scroll_offset_x = self.scroll_area.horizontalScrollBar().value()
            scroll_offset_y = self.scroll_area.verticalScrollBar().value()
            corrected_x = (event.pos().x() + scroll_offset_x) / self.scale_factor
            corrected_y = (event.pos().y() + scroll_offset_y) / self.scale_factor    
            self.lastDrawPoint = QPointF(corrected_x, corrected_y)
        else:
            self.dragging = True
            self.lastMousePosition = event.pos()
            self.setCursor(Qt.ClosedHandCursor) 

    def mouseMoveEvent(self, event):
        super().mouseMoveEvent(event)  # Call the parent class's method
        if self.isDrawing and self.drawingToolActive:
                currentTime = time.time()
                # First, get the current scroll positions
                scroll_offset_x = self.scroll_area.horizontalScrollBar().value()
                scroll_offset_y = self.scroll_area.verticalScrollBar().value()
                # Adjust the event position for the scroll offsets and the current zoom level
                # Since drawing is based on the original image, we scale the mouse position back
                # to "image space" by dividing by the zoom factor and adding the scroll offsets
                corrected_x = (event.pos().x() + scroll_offset_x)/ self.scale_factor
                corrected_y = (event.pos().y() + scroll_offset_y)/ self.scale_factor
                # Apply a manual correction factor to the Y-coordinate
                manual_correction_y = 50  # Adjust this value based on observation
                corrected_y -= manual_correction_y
                next_index = len(self.drawnPixels)
                self.drawnPixels.loc[next_index] = [corrected_x, corrected_y]                # Update lastDrawPoint for continuous drawing, using QPointF for more precision if needed
                if self.lastDrawPoint is None:
                    self.lastDrawPoint = QPointF(corrected_x, corrected_y)
                elif self.lastUpdateTime is None or currentTime - self.lastUpdateTime > 0.1:
                    currentPoint = QPointF(corrected_x, corrected_y)
                    # Perform the drawing on the drawablePixmap, which is scaled and/or scrolled appropriately
                    painter = QPainter(self.drawablePixmap)
                    pen = QPen(Qt.black, 3, Qt.SolidLine)  # Adjust pen properties as needed
                    painter.setPen(pen)
                    # Draw line from the last point to the current point
                    painter.drawLine(self.lastDrawPoint, currentPoint)
                    # Update lastDrawPoint for the next move event
                    self.lastDrawPoint = currentPoint
                    painter.end()  # End painting
                    # Update the image display to show the new drawing
                    self.updateImageDisplay()
        elif self.dragging and not self.drawingToolActive:
            self.setCursor(Qt.ClosedHandCursor)
            # calcute the distance the mouse has moved
            delta = event.pos() - self.lastMousePosition
            self.lastMousePosition = event.pos()
            # scroll the scroll area by the distance moved
            self.scroll_area.horizontalScrollBar().setValue(self.scroll_area.horizontalScrollBar().value() - delta.x())
            self.scroll_area.verticalScrollBar().setValue(self.scroll_area.verticalScrollBar().value() - delta.y())
            pass

    def mouseReleaseEvent(self, event):
        super().mouseReleaseEvent(event)  # Call the parent class's method
        if event.button() == Qt.LeftButton and self.drawingToolActive:
            self.isDrawing = False
        elif self.dragging:
            self.dragging = False
            self.setCursor(Qt.ArrowCursor) ##################

    def uploadCsv(self):
        # Let the user choose a CSV file
        csvFileName, _ = QFileDialog.getOpenFileName(self, 'Open CSV file', '', 'CSV files (*.csv)')
        if csvFileName:
            self.csvFilePath = csvFileName
            self.loadCsvData()
            self.backtrackToolAction.setEnabled(True)  # Enable backtrack tool after CSV is loaded

    def loadCsvData(self):
        # Read the specified CSV file and store the in_tissue == 1 positions
        if self.csvFilePath:
            df = pd.read_csv(self.csvFilePath)
             # Ask the user to input the scale factor
            text, ok = QtWidgets.QInputDialog.getText(
                self, 'Scale Factor', 'Enter the scale factor (e.g., 1.2345678901):'
            )
            if ok and text:
                try:
                    # Convert the input string to a float
                    scale_factor = float(text)
                    # Apply the scale factor and round the columns to the nearest integer
                    df['pxl_col_in_fullres'] = (df['pxl_col_in_fullres'] * scale_factor).round().astype(int)
                    df['pxl_row_in_fullres'] = (df['pxl_row_in_fullres'] * scale_factor).round().astype(int)
                except ValueError:
                    # Handle the case where the conversion fails
                    QtWidgets.QMessageBox.warning(self, 'Input Error', 'Invalid scale factor entered. Please enter a valid number.')
            self.tissue_positions_complete = df
            self.tissuePositions = df[df['in_tissue'] == 1][['pxl_col_in_fullres', 'pxl_row_in_fullres']].values


    def applyBacktrackTool(self):
        # Load tissue positions if not already loaded
        if self.tissuePositions is None:
            self.loadCsvData()

        if self.tissuePositions is not None and not self.originalPixmap.isNull():
            if not self.highlightLayer:
                self.highlightLayer = QPixmap(self.originalPixmap.size())
                self.highlightLayer.fill(Qt.transparent)
            
            painter = QPainter(self.highlightLayer)
            pen = QPen(Qt.black, 10)  # Adjust pen color and size as needed
            painter.setPen(pen)
            for position in self.tissuePositions:
                # Convert position to QPoint and draw
                painter.drawPoint(QPoint(position[0], position[1]))
            painter.end()

            # Update the display to show the highlights
            self.updateImageDisplay()

    
    #Implementing the mapping_script
    def mapping_script(self, dp):
        tissue_positions_df = self.tissue_positions_complete
        tissue_positions_df = tissue_positions_df[tissue_positions_df['in_tissue'] == 1]
        tree = KDTree(tissue_positions_df[['pxl_col_in_fullres', 'pxl_row_in_fullres']])
        closest_points_indices = tree.query(dp[['x', 'y']])[1]
        closest_barcodes = tissue_positions_df.iloc[closest_points_indices]['barcode']
        unique_barcodes = closest_barcodes.drop_duplicates()
        unique_barcodes.to_csv('barcodes.csv', index=False)
        print("Unique barcodes have been saved to barcodes.csv")

    #Implementing the extractData method
    def extractBarcodes(self):
        # Ensure that there are drawn pixels and the tissue positions have been loaded
        if self.tissuePositions is not None and not self.drawnPixels.empty:
            self.extractToolAction.setEnabled(True)
            # Call the mapping function with the DataFrame of drawn pixels
            barcodes = self.mapping_script(self.drawnPixels)
            # Here, handle the resulting barcodes, e.g., displaying them or logging
            print("Mapped Barcodes:", barcodes)
            # Optionally clear the drawnPixels DataFrame after processing
            self.drawnPixels = pd.DataFrame(columns=['x', 'y'])
        else:
            # Disable the extract tool if conditions are not met
            # self.extractToolAction.setEnabled(False)
            if self.tissuePositions is None:
                print("Tissue positions have not been loaded.")
            if self.drawnPixels.empty:
                print("No drawn pixels to process.")

    def clearAll(self):
        # Clear the drawable pixmap by filling it with a transparent color
        if self.drawablePixmap:
            self.drawablePixmap.fill(Qt.transparent)
            self.updateImageDisplay()
        # Optionally, clear other layers like highlightLayer if applicable
        if self.highlightLayer:
            self.highlightLayer.fill(Qt.transparent)
        # Clear stored pixel data
        self.drawnPixels = pd.DataFrame(columns=['x', 'y'])
        # Reset additional tool states if necessary
        self.drawToolAction.setChecked(False)
        self.drawingToolActive = False
        self.image_label.setCursor(Qt.ArrowCursor)

        print("All drawings and data have been cleared.")



    # setupTools method
    def setupTools(self):
        # draw tool setup
        self.drawToolAction = self.makeAction(self, './utils/paint.ico', 'Draw', 'Draw with pencil', self.activateDrawTool)
        self.drawToolAction.setCheckable(True)
        self.toolsmenu.addAction(self.drawToolAction)
        self.toolstoolbar.addAction(self.drawToolAction)

        # add more tools here
        self.uploadCsvAction = self.makeAction(self, './utils/upload.ico', 'Upload CSV', 'Upload tissue_position.csv', self.uploadCsv)
        self.backtrackToolAction = self.makeAction(self, './utils/backtrack.ico', 'Backtrack', 'Highlight tissue positions', self.applyBacktrackTool)
        self.backtrackToolAction.setEnabled(False)  # Disabled until a CSV is uploaded
        self.toolsmenu.addActions([self.uploadCsvAction, self.backtrackToolAction])
        self.toolstoolbar.addActions([self.uploadCsvAction, self.backtrackToolAction])

         # Clear tool setup
        self.clearToolAction = self.makeAction(self, './utils/clear.ico', 'Clear All', 'Clear all drawings and reset', self.clearAll)
        self.toolsmenu.addAction(self.clearToolAction)
        self.toolstoolbar.addAction(self.clearToolAction)

        # initialize drawing state
        self.isDrawing = False
        self.drawingToolActive = False
        self.drawToolAction.setEnabled(False)
    
    # activate draw tool
    def activateDrawTool(self):
        self.drawingToolActive = not self.drawingToolActive
        self.drawToolAction.setChecked(self.drawingToolActive)
        if self.drawingToolActive:
            self.image_label.setCursor(Qt.CrossCursor)
        else:
            # Reset cursor or other states when tool is deactivated
            self.image_label.setCursor(Qt.ArrowCursor)

    #defining the method to create the action for the menu options
    def makeAction(self, parent, icon, name, tip, method):  
        # creating an object of the QAction class to create the action for the menu options  
        the_action = QAction(QIcon(icon), name, parent)  
        # setting the status tip of the action using the setStatusTip() method  
        the_action.setStatusTip(tip)  
        # connecting the action to the method using the triggered.connect() method  
        the_action.triggered.connect(method)  
        # returning the action  
        return the_action

# main function  
if __name__ == '__main__':  

    # creating an object of the QApplication class  

    the_app = QApplication(sys.argv)  

    # creating an object of the Application class  
    imageViewerApp = QImageViewer()  

    # using the show() method to display the window  
    imageViewerApp.show()  

    # using the exit() function of the sys module to close the application  
    sys.exit(the_app.exec_())
