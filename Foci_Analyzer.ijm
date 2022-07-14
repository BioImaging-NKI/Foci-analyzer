#@ File[] (label = "Input files", style="File") files
#@ File (label = "Output folder", style = "directory") outputFolder

#@ Integer (label = "Nuclei channel", value = 1) nucleiChannel
#@ Integer (label = "Foci channel A", value = 2) fociChannelA
#@ Integer (label = "Foci channel B", value = 3) fociChannelB
#@ Boolean (label = "Also detect foci in channel B and perform colocalization?", persistence=true, value=true, description="Check this box to calculate colocalization between two foci channels.") detectFociChannelB

#@ Integer (label = "Image XY binning before analysis [0-n]", value = 1, min=1, description="Image binning will greatly shorten the analysis time") XYBinning

#@ String (value="Nuclei detection settings", visibility="MESSAGE") nuclei_message
#@ String (label = "Nuclei segmentation method", choices={"StarDist nuclei segmentation","Cellpose cytoplasm segmentation", "Classic nuclei segmentation"}, style="listBox") nucleiSegmentationChoice
#@ Integer (label = "Nuclei binning factor [0-n] (StarDist)", value = 2, min=1, description="A binning >1 for high-resolution input images can prevent 'oversegmentation', and will also somewhat speed up processing.") downsampleFactorStarDist
#@ Double (label = "Probablility threshold [0.0-1.0] (StarDist/Cellpose)", value = 0.5, min=0, max=1, description="lower values will accept more nuclei") probabilityThreshold

#@ Integer (label = "Remove nulei with diameter smaller than (units)", value = 4) minNucleusSize_setting
#@ Integer (label = "Remove nulei with diameter larger than (units)", value = 50) maxNucleusSize_setting
#@ Boolean (label = "Exclude nulei on image edges", value = true) excludeOnEdges
#@ String (label = "Manual nuclei removal", choices={"No thanks","Manually remove nuclei","load previously saved removals (from output folder)"}, value = "No thanks", description="Allow nuclei removal by clicking, or load previous removals (Masks should be in the output folder)") manualRemoveNuclei

#@ String (value="Foci detection settings", visibility="MESSAGE") foci_message1
#@ Boolean (label = "Enable foci parameters optimization mode?", value=true) optimizationMode
#@ String (label = "Foci size channel A", choices={"tiny","small","average","large","huge","other"}, style="listBox", value="average", description="Simplified foci detection parameter") fociSizeA
#@ String (label = "Foci size channel B", choices={"tiny","small","average","large","huge","other"}, style="listBox", value="average", description="Simplified foci detection parameter") fociSizeB
#@ String (label = "Foci detection method", choices={"Marker-controlled watershed","AreaMaxima local maximum detection"}, style="listBox") foci_method

#@ Double (label = "Foci intensity threshold bias channel A", value = 0.0, min=-2.5, max=2.5, style="scroll bar", stepSize=0.01) thresholdFactorA
#@ Double (label = "Foci intensity threshold bias channel B", value = 0.0, min=-2.5, max=2.5, style="scroll bar", stepSize=0.01) thresholdFactorB
#@ Double (label = "Minimum foci size (pixels/voxels)", value = 1) minFociSize
#@ Double (label = "Maximum foci size (pixels/voxels)", value = 9999) maxFociSize

#@ Boolean (label = "[3D only] Do foci detection on maximum intensity projection only?", value=false) maxProject
#@ Boolean (label = "[3D only] Use quasi-2D foci detection? (detect foci in every Z-slice)", value=false, description="For 3D stacks with large spacing between the z-slices") quasi2D
#@ Integer (label = "[3D only] Process a single z-slice only: slice nr (-1 for full stack)", value=-1, min=-1) singleSlice

#@ String (value="Colocalization options", visibility="MESSAGE") coloc_message
#@ Integer (label = "Minimum overlap of foci to colocalize (pixels/voxels)", min = 1, value = 1, description="If the overlap between foci in channel A and channel B is too small it will not be counted as colocalized") minOverlapSize

#@ Boolean (label = "Also save colocalization images", value=true) saveOverlayImage

#@ ColorRGB(label = "Cell number font color", value="orange") fontColor

#@ Boolean (label = "Debug mode (show intermediate images)", value=false) debugMode

version = 1.0;

/* Macro to quantify foci in nuclei/cells. More info will follow :-)
 * Works on 2D/3D images, including multiseries files, as long as all series have the same specs
 *   
 * ► Requires the following Fiji update sites:
 * - 3D ImageJ Suite
 * - CLIJ
 * - CLIJ2
 * - CLIJx-assistent
 * - CLIJx-assistent-extensions
 * - CSBDeep
 * - IJPB-plugins
 * - SCF MPI CBG
 * - StarDist
 * 
 * ► In order to run Cellpose segmentation you also need:
 * - A working Cellpose Python environment
 * - PTBIOP update site, with proper settings. See https://github.com/BIOP/ijl-utilities-wrappers/blob/master/README.md#cellpose
 * 
 * 
 *	Author: Bram van den Broek, The Netherlands Cancer Institute, b.vd.broek@nki.nl
 * 
 * 
 * 
 * Changelog
 * 
 * ---------
 * version 0.94, April 2022:
 * - Invert foci before MorpholibJ marker-controlled watershed -> get rid of holes in the foci.
 * - Convert foci to 32-bit before DoG filtering. Currently not converted back to 16 bit when determining foci thresholds.
 *   Both changes will substantially impact the measured foci from version 0.93!
 * 
 * Version 0.95, April 2022:
 * - Fixed overwriting of output files for multi-series images (e.g. .lif, .czi)
 * - Fixed mistakes in output metadata (foci size A and B were both A, 2 x 'threshold bias A' in the key)
 * - Added nucleus 2D area as output column
 * 
 * Version 0.96, April 2022:
 * - Added Cellpose as segmentation possibility (requires PTBIOP update site and a working Cellpose installation in Python)
 * 
 * Version 0.97, May 2022:
 * - Added possibility to bin the image in XY (all channels)
 * 
 * Version 0.98, June 2022:
 * - Fixed a bug: foci detection on max projections with 2 foci channels does not give an error any more
 * 
 * Version 0.99, June 2022:
 * - Updated the parameters saved in the metadata of the output images
 * 
 * * Version 1.0, July 2022:
 * - Added possibility to manually remove segmented nuclei by clicking. Masks are saved and can be loaded for re-analysis.
 */



/* TODO / IDEAS
 *  
 * Use ALT to pop up options menu (background subtraction, nuclear segmentation, foci size, etc.
 * and allow the user to determine the right settings! See MosaicExplorerJ macro.
 * Currently, pressing ALT during opening the image allows for making a selection.
 * 
 * Upscale and smooth ROIs (spline interpolate met Roi.setPolygonSplineAnchors(x, y)?)
 * 
 * Create log file with settings and/or save metadata in the output image
 * 
 * Add description="..." to the options (onmouseover)
 * 
 * Background intensity is now very crude (just the mean outside the nuclei). So: mask the background, set (percentile) threshold, measure mean/median. 
 * 
 */

//nuclei detection
nucleiBlurRadiusXY = 1;
nucleiBlurRadiusZ = 1;
nucleiMedian3DradiusXY = 2;
nucleiMedian3DradiusZ = 2;
minNucleusSize = 4;
maxNucleusSize = 40;
excludeOnEdges = excludeOnEdges;

//foci detection
var fociFilterSizeXY;	//foci radius for DifferenceOfGaussians
var fociFilterSizeZ;
var fociSizeXY;			//foci radius for detectMaximaBox
var fociSizeZ;

starDistTiles = 1;
thresholdMultiplier = 3;	//Default threshold multiplier
minFociSize = minFociSize;
maxFociSize = maxFociSize;

H_max_peak_flooding = 100;	//Percentage of the peak
thickOutlines = false;		//Width of nuclei outlines, 1 or 2 pixels
thickOutlinesThresholdSize = 1000;	//Above this image size nuclei outlines are thick
hideImages = false;

labelOpacity = 100;			//opacity of outlines
addNumbersOverlay = true;
labelFontSize = 12;

var processAllOtherImages = false;
var doneOptimizing;
var displayMode = "composite";
var activeChannels = "111";
var maxDisplaySetting;
var processTime = 0;
var threshold = 0;


//Create the azure LUT for the nuclei (TO DO: do this in a more elegant way than globals)
reds = newArray(256);
greens = newArray(256); 
blues = newArray(256);
create_azure_lut(reds, greens, blues);


saveSettings();

run("Set Measurements...", "area mean standard integrated median redirect=None decimal=3");
run("Conversions...", " ");
setOption("BlackBackground", true);
setForegroundColor(255, 255, 255);

roiManager("Reset");
run("Clear Results");
print("\\Clear");
close("\\Others");

run("Clear Results");
run("CLIJ2 Macro Extensions", "cl_device=");
Ext.CLIJ2_clear();


// Start of workflow
//------------------------------
setBatchMode(true);

run("Close All");
nrOfImages = files.length;
for (f = 0; f < nrOfImages; f++) {
	print("\nProcessing file "+f+1+"/"+nrOfImages+": "+files[f]);
	processFile(f, files[f], outputFolder);
}
print("\nFinished processing "+nrOfImages+" files in "+processTime*60+" seconds ("+d2s(processTime,1)+" minutes).");
print("\\Update3:Finished processing "+nrOfImages+" files in "+processTime*60+" seconds ("+d2s(processTime,1)+" minutes).");

restoreSettings();


function processFile(current_image_nr, file, outputFolder) {
	run("Close All");

	startTimeSeries = getTime();
	print("\\Update1:Processing file "+current_image_nr+1+"/"+nrOfImages+": " + file);
	print("\\Update2:Average speed: "+d2s((current_image_nr)/processTime,1)+" images per minute.");
	time_to_run = (nrOfImages-(current_image_nr)) * processTime/(current_image_nr);
	if(time_to_run<5) print("\\Update3:Projected run time: "+d2s(time_to_run*60,0)+" seconds ("+d2s(time_to_run,1)+" minutes).");
	else if(time_to_run<60) print("\\Update3:Projected run time: "+d2s(time_to_run,1)+" minutes. You'd better get some coffee.");
	else if(time_to_run<480) print("\\Update3:Projected run time: "+d2s(time_to_run,1)+" minutes ("+d2s(time_to_run/60,1)+" hours). You'd better go and do something useful.");
	else if(time_to_run<1440) print("\\Update3:Projected run time: "+d2s(time_to_run,1)+" minutes. ("+d2s(time_to_run/60,1)+" hours). You'd better come back tomorrow.");
	else if(time_to_run>1440) print("\\Update3:Projected run time: "+d2s(time_to_run,1)+" minutes. This is never going to work. Give it up!");
	print("\\Update4:-------------------------------------------------------------------------");

	run("Bio-Formats Macro Extensions");	//Necessary to do this here, because you can only activate one Macro Extension at the time
	Ext.setId(file);
	Ext.getSeriesCount(nr_series);

	run("CLIJ2 Macro Extensions", "cl_device="); //Necessary to do this here, because you can only activate one Macro Extension at the time
	Ext.CLIJ2_clear();

	filename = File.getName(file);
	fileExtension = substring(filename, lastIndexOf(filename, "."), lengthOf(filename));
	if(endsWith(fileExtension, "tif") || endsWith(fileExtension, "jpg") || endsWith(fileExtension, "png")) {	//Use standard opener
		open(file);
		if(XYBinning > 1) run("Bin...", "x="+XYBinning+" y="+XYBinning+" z=1 bin=Average");
		process_current_series(file, true);
	}
	else {	//Use Bio-Formats
		for(s = 0; s < nr_series; s++) {
			run("Close All");
			run("Bio-Formats Importer", "open=["+file+"] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT series_"+s+1);
			seriesName = getTitle();
			seriesName = replace(seriesName,"\\/","-");	//replace slashes by dashes in the seriesName
			print(seriesName);
	//		outputPath = output + File.separator + substring(seriesName)

			if(XYBinning > 1) run("Bin...", "x="+XYBinning+" y="+XYBinning+" z=1 bin=Average");
			process_current_series(seriesName, false);
		}
	}
}



function process_current_series(image, nameHasExtension) {
	if(singleSlice != -1) run("Duplicate...", "duplicate slices="+singleSlice);
	Stack.setDisplayMode("grayscale");
	Stack.setChannel(nucleiChannel);
	original = getTitle();

	if(nameHasExtension) imageName = File.getNameWithoutExtension(image);
	else imageName = image;
	
	//Initialize image and table
	getDimensions(gwidth, gheight, gchannels, gslices, gframes); // global variables
	getVoxelSize(pixelWidth, pixelHeight, pixelDepth, unit);
	if(gwidth > thickOutlinesThresholdSize && gheight > thickOutlinesThresholdSize) thickOutlines = true;
	
	if(gslices > 1) anisotropyFactor = pixelWidth / pixelDepth;
	else anisotropyFactor = 0;
	bits = bitDepth();
	if(bits == 8) run("16-bit");	//Convert to 16-bit, because 8-bit restricts the foci labelmap to 255 foci

	Stack.setSlice(gslices/2);
	setBatchMode("show");
	run("Enhance Contrast", "saturated=0.35");

	resultTable = "Foci results";
	if(isOpen(resultTable)) close(resultTable);
	Table.create(resultTable);
	Table.setLocationAndSize(0, 0, 1000, 500);

	if(gslices > 1 && maxProject == true) {
		selectWindow(original);
		run("Z Project...", "projection=[Max Intensity]");
		original = getTitle();
		gslices = 1;
	}

	//Segment and label nuclei
	if(nucleiSegmentationChoice == "Classic nuclei segmentation") {
		nuclei_info = segmentNuclei(original, nucleiChannel, nucleiBlurRadiusXY, nucleiBlurRadiusZ, nucleiMedian3DradiusXY, nucleiMedian3DradiusZ);
	}
	else if (nucleiSegmentationChoice == "StarDist nuclei segmentation") nuclei_info = segmentNucleiStarDist(original, nucleiChannel, probabilityThreshold, pixelWidth, unit, resultTable);
	else if (nucleiSegmentationChoice == "Cellpose cytoplasm segmentation") nuclei_info = segmentCytoplasmCellpose(original, nucleiChannel, probabilityThreshold, pixelWidth, unit, resultTable);
	labelmap_nuclei = nuclei_info[0];
	nuclei_outlines = nuclei_info[1];
	nrNuclei = nuclei_info[2];

	if(manualRemoveNuclei != "No thanks" && nrNuclei > 0) {
		nuclei_info = manually_remove_labels(labelmap_nuclei, nuclei_outlines, nrNuclei, original, imageName);
		labelmap_nuclei = nuclei_info[0];
		nuclei_outlines = nuclei_info[1];
		nrNuclei = nuclei_info[2];
		labelmap_nuclei = "Labelmap_nuclei_edited";
	}

	if(nrNuclei > 0) {
		//Write ID and area to foci results table
		run("Clear Results");
		Ext.CLIJ2_statisticsOfLabelledPixels(labelmap_nuclei, labelmap_nuclei);
		nucleus_id_ = Table.getColumn("IDENTIFIER", "Results");
		nucleus_area_ = Table.getColumn("PIXEL_COUNT", "Results");
		nucleus_area_ = multiplyArraywithScalar(nucleus_area_, pixelWidth / Math.sqr(downsampleFactorStarDist));
		Table.setColumn("Nucleus ID", nucleus_id_, resultTable);
		Table.setColumn("Nucleus area 2D ("+unit+"^2)", nucleus_area_, resultTable);
	
		//Create a 3D version of the 2D nuclei labelmap
		labelmap_nuclei_3D = createLabelmapNuclei3D(labelmap_nuclei);
	
		//Foci filtering and detection - in a loop to enable parameter optimization 
		firstTimeProcessing = true;
	//	activeChannels = "111";
		zoom = 1;
		do {
			//Foci channel A
//			filtered_fociA = filterFoci(original, fociChannelA, fociSizeA, anisotropyFactor, ROFtheta, firstTimeProcessing);
			labelmap_fociA = detectFoci(original, fociChannelA, fociSizeA, anisotropyFactor, firstTimeProcessing, labelmap_nuclei, labelmap_nuclei_3D, thresholdFactorA);
			if(firstTimeProcessing == false) call("ij.gui.ImageWindow.setNextLocation", x_image, y_image);
			mergeOriginalAndDetection(original, nrNuclei, nuclei_outlines, "Mask_foci_ch"+fociChannelA, "foci_spots_ch"+fociChannelA, fociChannelA);
			overlay_image = getTitle();
			if(firstTimeProcessing == false && detectFociChannelB == false) close("Processing...");

			//Foci channel B
			if(detectFociChannelB == true) {
//				filtered_fociB = filterFoci(original, fociChannelB, fociSizeB, anisotropyFactor, ROFtheta, firstTimeProcessing);
				labelmap_fociB = detectFoci(original, fociChannelB, fociSizeB, anisotropyFactor, firstTimeProcessing, labelmap_nuclei, labelmap_nuclei_3D, thresholdFactorB);
				mergeOriginalAndDetection(original, nrNuclei, nuclei_outlines, "Mask_foci_ch"+fociChannelB, "foci_spots_ch"+fociChannelB, fociChannelB);
				if(firstTimeProcessing == false) close("Processing...");
				Overlay.copy();	//Preserve nuclei outline overlays
				overlay_image = "Foci_overlay_ch"+fociChannelA+"_and_ch"+fociChannelB;
				run("Concatenate...", "  title="+overlay_image+" image1=Foci_overlay_ch"+fociChannelA+" image2=Foci_overlay_ch"+fociChannelB+" image3=[-- None --]");
				run("Stack to Hyperstack...", "order=xyczt(default) channels=4 slices="+gslices+" frames=2 display=Composite");
				Overlay.paste();
				Stack.setFrame(2);
				run("Label...", "format=Text starting=0 interval=0 x=20 y=10 font=18 text=[Channel "+fociChannelB+"] range=[2-2] use");
				Stack.setFrame(1);
				run("Label...", "format=Text starting=0 interval=0 x=20 y=10 font=18 text=[Channel "+fociChannelA+"] range=[1-1] use");
			}
			//Set display settings and location
			if(firstTimeProcessing == false) {
				Stack.setChannel(1);
				setMinAndMax(minDisplayNuclei, maxDisplayNuclei);
				Stack.setChannel(2);
				setMinAndMax(minDisplayFoci, maxDisplayFoci);
//				Stack.setDisplayMode(displayMode);
//				Stack.setActiveChannels(activeChannels);
			}
			Stack.setDisplayMode(displayMode);
			Stack.setActiveChannels(activeChannels);
			if(firstTimeProcessing == false) Stack.setPosition(currentChannel, currentSlice, currentFrame);
			else if(gslices > 1) Stack.setSlice(gslices/2);
			setBatchMode("show");

			if(firstTimeProcessing == false) run("Set... ", "zoom="+zoom*100+" x="+displayX + displayWidth/2+" y="+displayY + displayHeight/2);
			run("Channels Tool...");

			if(firstTimeProcessing == true) getLocationAndSize(x_image, y_image, imageWidth, imageHeight);

			//Parameter optimization dialog
			if(optimizationMode == true && processAllOtherImages == false) {

				Dialog.createNonBlocking("Optimize settings for foci detection");
				Dialog.addChoice("Foci size channel "+fociChannelA, newArray("tiny","small","average","large","huge","other"), fociSizeA);
				if(detectFociChannelB) Dialog.addChoice("Foci size channel "+fociChannelB, newArray("tiny","small","average","large","huge","other"), fociSizeB);
				Dialog.addChoice("Detection method", newArray("Marker-controlled watershed","AreaMaxima local maximum detection"), foci_method);
				Dialog.addSlider("Threshold bias channel "+fociChannelA+" (higher is more strict)", -2.5, 2.5, thresholdFactorA);
				if(detectFociChannelB) Dialog.addSlider("Threshold bias channel "+fociChannelB+" (higher is more strict)", -2.5, 2.5, thresholdFactorB);
				if(gslices > 1) {
					Dialog.addNumber("Minimum foci size", minFociSize, 1, 4, "voxels");
					Dialog.addNumber("Maximum foci size", maxFociSize, 1, 4, "voxels");
				}
				else {
					Dialog.addNumber("Minimum foci size", minFociSize, 1, 4, "pixels");
					Dialog.addNumber("Maximum foci size", maxFociSize, 1, 4, "pixels");
				}
				if(detectFociChannelB) {
					if(gslices > 1) Dialog.addNumber("Minimum overlap of foci", minOverlapSize, 1, 4, "voxels");
					else Dialog.addNumber("Minimum overlap of foci", minOverlapSize, 1, 4, "pixels");
				}
				Dialog.addCheckbox("Done optimizing for this image", false);
				Dialog.addCheckbox("Process all other images with these settings", false);

				//determine dialog location
				if(x_image + imageWidth + 500 < screenWidth) Dialog.setLocation(x_image+imageWidth-15, y_image);
				else if(x_image > 500) Dialog.setLocation(x_image-530, y_image);
//				else if(y_image + gheight + 300 < screenHeight) Dialog.setLocation(x_image, y_image+gheight+100);
				else Dialog.setLocation(x_image, y_image);
				Dialog.show();

				//Get dialog entries
				fociSizeA = Dialog.getChoice();
				if(detectFociChannelB == true) fociSizeB = Dialog.getChoice();
				foci_method = Dialog.getChoice();
				thresholdFactorA = Dialog.getNumber();
				if(detectFociChannelB) thresholdFactorB = Dialog.getNumber();
				minFociSize = Dialog.getNumber();
				maxFociSize = Dialog.getNumber();
				if(detectFociChannelB) minOverlapSize = Dialog.getNumber();
				doneOptimizing = Dialog.getCheckbox();
				processAllOtherImages = Dialog.getCheckbox();
				if (processAllOtherImages == true) doneOptimizing = true;

				//get image display properties
				getLocationAndSize(x_image, y_image, imageWidth, imageHeight);
				Stack.getActiveChannels(activeChannels);
				Stack.getDisplayMode(displayMode);
				Stack.setChannel(1);
				getMinAndMax(minDisplayNuclei, maxDisplayNuclei);
				Stack.setChannel(2);
				getMinAndMax(minDisplayFoci, maxDisplayFoci);
				Stack.getPosition(currentChannel, currentSlice, currentFrame);
				zoom = getZoom();
				getDisplayedArea(displayX, displayY, displayWidth, displayHeight);

				if(doneOptimizing == false) {
					close("foci_ch"+fociChannelA);
					//close("MAX_Foci_ch"+channel+"_filtered_and_masked");
					close("Labelmap_detected_foci_filtered_ch"+fociChannelA);
					close("Mask_foci_ch"+fociChannelA);
//					close("foci_RAW"+fociChannelA);
					if(detectFociChannelB == false) {
						selectWindow("Foci_overlay_ch"+fociChannelA);
						rename("Processing...");
					}
					else if(detectFociChannelB == true) {
						close("foci_ch"+fociChannelB);
						//close("MAX_Foci_ch"+channel+"_filtered_and_masked");
						close("Labelmap_detected_foci_filtered_ch"+fociChannelB);
						close("Mask_foci_ch"+fociChannelB);
//						close("foci_RAW_ch"+fociChannelB);
						selectWindow("Foci_overlay_ch"+fociChannelA+"_and_ch"+fociChannelB);
						rename("Processing...");
					}
				}

				// Create table with settings
				Table.create("Last foci detection settings");
				Table.set("Parameter", Table.size, "Foci size channel "+fociChannelA);
				Table.set("Value", Table.size-1, fociSizeA);
				Table.set("Parameter", Table.size, "Foci size channel "+fociChannelB);
				Table.set("Value", Table.size-1, fociSizeB);
//					Table.set("Value", 0, d2s(fociFilterSizeXY, 2));
//					Table.set("Parameter name", 0, "Foci filter size xy");
//					Table.set("Value", 0, d2s(fociFilterSizeXY, 2));
//					Table.set("Parameter name", 1, "Foci filter size z");
//					Table.set("Value", 1, d2s(fociFilterSizeZ, 2));
//					Table.set("Parameter name", 2, "Foci detection size xy");
//					Table.set("Value", 2, d2s(fociSizeXY, 2));
//					Table.set("Parameter name", 3, "Foci detection size z");
//					Table.set("Value", 3, d2s(fociSizeZ, 2));
				Table.set("Parameter", Table.size, "Detection method");
				Table.set("Value", Table.size-1, foci_method);
				Table.set("Parameter", Table.size, "Threshold channel A");
				Table.set("Value", Table.size-1, d2s(thresholdFactorA, 2));		
				Table.set("Parameter", Table.size, "Threshold channel B");
				Table.set("Value", Table.size-1, d2s(thresholdFactorB, 2));
				Table.set("Parameter", Table.size, "Minimum foci size");
				Table.set("Value", Table.size-1, d2s(minFociSize, 0));		
				Table.set("Parameter", Table.size, "Maximum foci size");
				Table.set("Value", Table.size-1, d2s(maxFociSize, 0));
				Table.update();
			}
//		Stack.setChannel(2);
		firstTimeProcessing = false;
		} while(optimizationMode == true && doneOptimizing == false);

		//Get nuclei positions
		selectWindow("Results");
		boundingBox_X = Table.getColumn("BOUNDING_BOX_X");
		boundingBox_Y = Table.getColumn("BOUNDING_BOX_Y");
		boundingBox_width = Table.getColumn("BOUNDING_BOX_WIDTH");
		boundingBox_height = Table.getColumn("BOUNDING_BOX_HEIGHT");
		Array.getStatistics(boundingBox_width, minWidth, maxWidth);
		Array.getStatistics(boundingBox_height, minHeight, maxHeight);
		print("Maximum nucleus bounding box: "+maxWidth+", "+maxHeight);		

		//Measure the foci 
		measureFoci(original, fociChannelA, nrNuclei, labelmap_nuclei_3D, labelmap_fociA, boundingBox_X, boundingBox_Y, maxWidth, maxHeight);
		if(detectFociChannelB) measureFoci(original, fociChannelB, nrNuclei, labelmap_nuclei_3D, labelmap_fociB, boundingBox_X, boundingBox_Y, maxWidth, maxHeight);

//		print("\n\nGPU Memory after channel "+fociChannelB);
//		Ext.CLIJ2_reportMemory();
	
		//Compute A-B foci colocalization
		if(detectFociChannelB) {
			computeOverlap("Mask_foci_ch" + fociChannelA, "Mask_foci_ch" + fociChannelB, nrNuclei, labelmap_nuclei_3D, boundingBox_X, boundingBox_Y, maxWidth, maxHeight);	//These images are still open in RAM/GPU
//			print("\n\nGPU Memory after computing overlap");
//			Ext.CLIJ2_reportMemory();
		}

		//Measure intensity in non-foci channels
		measure_nuclear_intensities(original, nrNuclei, labelmap_nuclei_3D, gchannels, fociChannelA, fociChannelB, resultTable);

		//Save data

		//Create parameter list
		getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
		MonthNames = newArray("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec");
     	DayNames = newArray("Sun", "Mon","Tue","Wed","Thu","Fri","Sat");
		List.set("00. =Foci Analyzer settings", "");
		List.set("00. Date", " " + DayNames[dayOfWeek] + " " + dayOfMonth + " " + MonthNames[month] + " " + year);
		List.set("00. Time", " " + hour +":"+IJ.pad(minute,2)+":"+IJ.pad(second,2));
		List.set("01. Version", version);
		List.set("02. Nuclei channel", nucleiChannel);
		List.set("03. Foci channel A", fociChannelA);
		List.set("04. Foci channel B", fociChannelB);
		List.set("05. Also detect foci channel B?", detectFociChannelB);
		List.set("06. Nuclei segmentation method", nucleiSegmentationChoice);
		List.set("07. Nuclei binning factor [0-n] (StarDist)", downsampleFactorStarDist);
		List.set("08. Nuclei probablility threshold [0.0-1.0] (StarDist)", probabilityThreshold);
		List.set("09. Remove nulei with diameter smaller than (units)", minNucleusSize_setting);
		List.set("10. Remove nulei with diameter larger than (units)", maxNucleusSize_setting);
		List.set("11. Exclude nulei on image edges", excludeOnEdges);
		List.set("12. Enable foci parameters optimization mode?", optimizationMode);
		List.set("13. Foci size channel A", fociSizeA);
		List.set("14. Foci size channel B", fociSizeB);
		List.set("15. Foci detection method", foci_method);
		List.set("16. Foci intensity threshold bias channel A", thresholdFactorA);
		List.set("17. Foci intensity threshold bias channel B", thresholdFactorB);
		List.set("18. Minimum foci size", minFociSize);
		List.set("19. Maximum foci size", maxFociSize);
		List.set("20. Do foci detection on maximum intensity projection only?", maxProject);
		List.set("21. Use quasi-2D foci detection? (detect foci in every Z-slice)", quasi2D);
		List.set("22. Process a single z-slice only: slice nr (-1 for full stack)", singleSlice);
		List.set("23. Minimum overlap of foci to colocalize (pixels/voxels)", minOverlapSize);
		list = List.getList();

		Table.save(outputFolder + File.separator + imageName + "__Foci_results.tsv");
		selectWindow(overlay_image);
		setMetadata("info", list);
		saveAs("tif", outputFolder + File.separator + imageName + "__" + overlay_image);
		if(saveOverlayImage == true) {
			if(detectFociChannelB == true) {
				selectWindow("Foci overlap map");
				setMetadata("info", list);
				saveAs("png", outputFolder + File.separator + imageName + "__foci_overlap_map");
			}
		}
		run("Clear Results");
		Ext.CLIJ2_release("nuclei");	
		Ext.CLIJ2_release(labelmap_nuclei);
		Ext.CLIJ2_release("foci_ch"+fociChannelA);
		if(detectFociChannelB == true) Ext.CLIJ2_release("foci_ch"+fociChannelB);
		Ext.CLIJ2_release("Labelmap_detected_foci_ch"+fociChannelA);
		if(detectFociChannelB == true) Ext.CLIJ2_release("Labelmap_detected_foci_ch"+fociChannelB);
//		Ext.CLIJ2_reportMemory();
	}
	else print("WARNING: 0 nuclei detected in "+ image +" !");
	
	//Ext.CLIJx_organiseWindows(Number startX, Number startY, Number tilesX, Number tilesY, Number tileWidth, Number tileHeight);
	
	cleanup();
	//Ext.CLIJ2_reportMemory();
	
	endTime = getTime();
	processTime = processTime+(endTime-startTimeSeries)/60000;
	//------------------------------
	// End of workflow
}



function segmentNuclei(image, channel, nucleiBlurRadiusXY, nucleiBlurRadiusZ, nucleiMedian3DradiusXY, nucleiMedian3DradiusZ) {
	selectWindow(image);
	getDimensions(width, height, channels, slices, frames);
	minNucleusSize = PI*pow((minNucleusSize_setting / pixelWidth / 2),2);	//Calculate the nucleus area as if it were a circle
	maxNucleusSize = PI*pow((maxNucleusSize_setting / pixelWidth / 2),2);
	
	nuclei = "nuclei";
	run("Duplicate...", "title=nuclei duplicate channels=&channel");
	run("32-bit");
	showStatus("Segmenting nuclei...");
	Ext.CLIJ2_push(nuclei);

	Ext.CLIJ2_gaussianBlur3D(nuclei, nuclei_filtered1, nucleiBlurRadiusXY, nucleiBlurRadiusXY, nucleiBlurRadiusZ);
	//run("Median 3D...", "x=" + nucleiMedian3Dradius + " y=" + nucleiMedian3Dradius + " z=" + nucleiMedian3Dradius);

//	CLIJ2_Median3D gives problems, at least on my GPU
//	Ext.CLIJ2_median3DSphere(nuclei_filtered1, nuclei_filtered2, nucleiMedian3DradiusXY/2, nucleiMedian3DradiusXY/2, nucleiMedian3DradiusZ/2);

	Ext.CLIJ2_maximumZProjection(nuclei_filtered1, nuclei_filtered_MAX);

	Ext.CLIJ2_automaticThreshold(nuclei_filtered_MAX, thresholded, "Huang");
	Ext.CLIJ2_pullBinary(thresholded);
	run("Watershed");
	run("Properties...", "unit=&unit pixel_width=&pixelWidth pixel_height=&pixelHeight voxel_depth=1.0000");
	if(excludeOnEdges) run("Analyze Particles...", "size=" + minNucleusSize + "-" + maxNucleusSize + " circularity=0.25-1.00 show=[Count Masks] exclude include add");
	else run("Analyze Particles...", "size=" + minNucleusSize + "-" + maxNucleusSize + " circularity=0.25-1.00 show=[Count Masks] include add");
	rename("nucleiCountMask");
	Ext.CLIJ2_push("nucleiCountMask");

	close("nuclei_filtered_MAX_blur");
	nrNuclei = roiManager("count");

	selectWindow(image);
	roiManager("Show All without labels");
	run("ROI Manager to LabelMap(2D)");
	rename("Labelmap_nuclei");
	return newArray("Labelmap_nuclei", nrNuclei);
}


function segmentNucleiStarDist (image, channel, probabilityThreshold, pixelWidth, unit, resultTable) {
	selectWindow(image);
	getDimensions(width, height, channels, slices, frames);
	minNucleusSize = PI*Math.sqr((minNucleusSize_setting / pixelWidth / 2));	//Calculate the nucleus area as if it were a circle
	maxNucleusSize = PI*Math.sqr((maxNucleusSize_setting / pixelWidth / 2));
	
	run("Duplicate...", "title=nuclei duplicate channels=" +channel);	//Get nucleus channel
	starDist_input_image = getTitle();
	
	//Downsample and Z-project.**** TO DO: Perform this on the GPU (but uses more GPU RAM)
	if(downsampleFactorStarDist > 1) {
		run("Scale...", "x="+1/downsampleFactorStarDist+" y="+1/downsampleFactorStarDist+" interpolation=None average process create");
		rename("nuclei_downscaled");
		starDist_input_image = getTitle();
	}
	if(slices>1) {
		run("Z Project...", "projection=[Max Intensity]");
		rename("nuclei_Zprojected");
		starDist_input_image = getTitle();
		if(downsampleFactorStarDist > 1) close("nuclei_downscaled");
	}

	//Run StarDist
	run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':'"+starDist_input_image+"', 'modelChoice':'Versatile (fluorescent nuclei)', 'normalizeInput':'true', 'percentileBottom':'1.0', 'percentileTop':'99.60000000000001', 'probThresh':'"+probabilityThreshold+"', 'nmsThresh':'0.3', 'outputType':'ROI Manager', 'nTiles':'"+starDistTiles+"', 'excludeBoundary':'2', 'roiPosition':'Stack', 'verbose':'false', 'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");
	run("ROI Manager to LabelMap(2D)");
	rename("Label Image");
	selectWindow("Label Image");
	run("Clear Results");

	//Scale up again
	if(downsampleFactorStarDist > 1) run("Scale...", "x="+downsampleFactorStarDist+" y="+downsampleFactorStarDist+" interpolation=None average process create");
	rename("Labelmap_nuclei_unfiltered");
	getDimensions(width, height, channels, slices, frames);
	if(width != gwidth || height != gheight) run("Canvas Size...", "width="+gwidth+" height="+gheight+" position=Center zero");	// Make sure that the size of the upscaled labelmap is correct
	labelmap_nuclei = "Labelmap_nuclei_unfiltered";

	//Close unused images
	close("Label Image");
	if(slices>1) close("nuclei_Zprojected");

	//Exclude nuclei on edges and count nr of nuclei
	Ext.CLIJ2_push(labelmap_nuclei);
	close(labelmap_nuclei);
	//Ext.CLIJ2_excludeLabelsWithValuesOutOfRange(area, labelmap_nuclei, labelmap_nuclei_filtered, minNucleusSize, maxNucleusSize);
	if(excludeOnEdges) {
		Ext.CLIJ2_excludeLabelsOnEdges(labelmap_nuclei, labelmap_nuclei_edges_excluded);
		Ext.CLIJ2_release(labelmap_nuclei);
	}
	else labelmap_nuclei_edges_excluded = labelmap_nuclei;
	Ext.CLIJ2_excludeLabelsOutsideSizeRange(labelmap_nuclei_edges_excluded, labelmap_nuclei_filtered, minNucleusSize, maxNucleusSize);
	Ext.CLIJ2_release(labelmap_nuclei_edges_excluded);
	labelmap_nuclei_final = "Labelmap_nuclei";
	Ext.CLIJ2_closeIndexGapsInLabelMap(labelmap_nuclei_filtered, labelmap_nuclei_final);
	Ext.CLIJ2_release(labelmap_nuclei_filtered);
	Ext.CLIJ2_getMaximumOfAllPixels(labelmap_nuclei_final, nrNuclei);	//Count the number of nuclei

	//Create nuclei outlines from nuclei labelmap
	Ext.CLIJ2_detectLabelEdges(labelmap_nuclei_final, labelmap_edges);
	if(thickOutlines == false) {
		Ext.CLIJ2_mask(labelmap_nuclei_final, labelmap_edges, labelmap_outlines);
		Ext.CLIJ2_release(labelmap_edges);
	}
	else labelmap_outlines = labelmap_edges;
	Ext.CLIJ2_pullBinary(labelmap_outlines);
	Ext.CLIJ2_release(labelmap_outlines);
	rename("nuclei_outlines");
	run("Cyan");
	
	nuclei_outlines = "nuclei_outlines";
	labelmap_nuclei_final = "Labelmap_nuclei";

	return newArray(labelmap_nuclei_final, nuclei_outlines, nrNuclei);
}



function segmentCytoplasmCellpose (image, channel, probabilityThreshold, pixelWidth, unit, resultTable) {
	selectWindow(image);
	getDimensions(width, height, channels, slices, frames);
	minNucleusSize = PI*Math.sqr((minNucleusSize_setting / pixelWidth / 2));	//Calculate the nucleus area as if it were a circle
	maxNucleusSize = PI*Math.sqr((maxNucleusSize_setting / pixelWidth / 2));

	run("Duplicate...", "title=nuclei duplicate channels=" +channel);	//Get nucleus channel
	setBatchMode("show");

	setBatchMode(false);
	run("Cellpose Advanced", "diameter=0 cellproba_threshold=0 flow_threshold="+probabilityThreshold+" anisotropy=1.0 diam_threshold=12.0 model=cyto nuclei_channel=0 cyto_channel=1 dimensionmode=2D stitch_threshold=-1.0 omni=false cluster=false additionnal_flags=-");
	labelmap_nuclei = getTitle();
//	selectWindow(image);
	selectWindow(labelmap_nuclei);
	setBatchMode("hide");
	selectWindow("nuclei");
	setBatchMode("hide");

	
	//Exclude nuclei on edges and count nr of nuclei
	Ext.CLIJ2_push(labelmap_nuclei);
	close(labelmap_nuclei);
	//Ext.CLIJ2_excludeLabelsWithValuesOutOfRange(area, labelmap_nuclei, labelmap_nuclei_filtered, minNucleusSize, maxNucleusSize);
	if(excludeOnEdges) {
		Ext.CLIJ2_excludeLabelsOnEdges(labelmap_nuclei, labelmap_nuclei_edges_excluded);
		Ext.CLIJ2_release(labelmap_nuclei);
	}
	else labelmap_nuclei_edges_excluded = labelmap_nuclei;
	Ext.CLIJ2_excludeLabelsOutsideSizeRange(labelmap_nuclei_edges_excluded, labelmap_nuclei_filtered, minNucleusSize, maxNucleusSize);
	Ext.CLIJ2_release(labelmap_nuclei_edges_excluded);
	labelmap_nuclei_final = "Labelmap_nuclei";
	Ext.CLIJ2_closeIndexGapsInLabelMap(labelmap_nuclei_filtered, labelmap_nuclei_final);
	Ext.CLIJ2_release(labelmap_nuclei_filtered);
	Ext.CLIJ2_getMaximumOfAllPixels(labelmap_nuclei_final, nrNuclei);	//Count the number of nuclei

	//Create nuclei outlines from nuclei labelmap
	Ext.CLIJ2_detectLabelEdges(labelmap_nuclei_final, labelmap_edges);
	if(thickOutlines == false) {
		Ext.CLIJ2_mask(labelmap_nuclei_final, labelmap_edges, labelmap_outlines);
		Ext.CLIJ2_release(labelmap_edges);
	}
	else labelmap_outlines = labelmap_edges;
	Ext.CLIJ2_pullBinary(labelmap_outlines);
	Ext.CLIJ2_release(labelmap_outlines);
	rename("nuclei_outlines");
	run("Cyan");
	
	nuclei_outlines = "nuclei_outlines";
	labelmap_nuclei_final = "Labelmap_nuclei";

	return newArray(labelmap_nuclei_final, nuclei_outlines, nrNuclei);
}


//Manually remove labels by creating point selections. Masks are saved as PNG files (compressed) for later use.
function manually_remove_labels(labelmap, label_edges, nrNuclei, windowName, imageName) {
	selectWindow(windowName);
	Stack.setDisplayMode("composite");
	Stack.setChannel(nucleiChannel);
    setLut(reds, greens, blues);
	Stack.setChannel(fociChannelA);
	run("Green");
   	if(detectFociChannelB) {
   		Stack.setChannel(fociChannelB);
   		run("Red");
   	}
   	run("Add Image...", "image="+label_edges+" x=0 y=0 opacity=100 zero");

	if(manualRemoveNuclei == "load previously saved removals (from output folder)" && File.exists(outputFolder + File.separator + imageName + "__removalMask.png")) {
		open(outputFolder + File.separator + imageName + "__removalMask.png");
		if(is("Inverting LUT")) run("Invert LUT");
		run("Points from Mask");
	}
	else if(manualRemoveNuclei == "Manually remove nuclei") {
		setTool("multipoint");
		run("Select None");
		run("Point Tool...", "type=Dot color=White size=[Extra Large] show counter=0");
		waitForUser("Select nuclei to remove, then press OK");
		if(selectionType == 10) {	//Only create mask if there is a selection
			run("Create Mask");
			saveAs("png", outputFolder + File.separator + imageName+"__removalMask");
		}
		else if(File.exists(outputFolder + File.separator + imageName + "__removalMask.png")) File.delete(outputFolder + File.separator + imageName + "__removalMask.png");	//Remove Mask
	}
	selectWindow(windowName);
	run("Select None");

	close(imageName + "__removalMask");

	Ext.CLIJ2_getMaximumOfAllPixels(labelmap, nrLabels);
	Ext.CLIJ2_pull(labelmap);
	selectWindow(labelmap);
	run("Restore Selection");
	run("Clear Results");
	run("Set Measurements...", "mean redirect=None decimal=3");
	run("Measure");
	run("Select None");
	setTool("rectangle");

	nuclei_to_be_removed = Table.getColumn("Mean", "Results");	//List of nuclei to be removed
	flaglist = "flaglist";
	newImage(flaglist, "32-bit black", nrLabels+1, 1, 1);	//+1 because background = 0
	for (i = 0; i < nuclei_to_be_removed.length; i++) setPixel(nuclei_to_be_removed[i], 0, 1);
	Ext.CLIJ2_push(flaglist);
	setBatchMode("show");
	close(flaglist);
	labelmap_edited = "Labelmap_nuclei_edited";
	Ext.CLIJ2_excludeLabels(flaglist, labelmap, labelmap_edited);
	Ext.CLIJ2_release(labelmap);
	labelmap = "Labelmap_nuclei_edited";

	Ext.CLIJ2_getMaximumOfAllPixels(labelmap_edited, nrNuclei);

	label_edges_edited = "label_edges_edited";
	Ext.CLIJ2_detectLabelEdges(labelmap_edited, label_edges_edited);
	Ext.CLIJ2_pull(label_edges_edited);
	run("Cyan");
	resetMinAndMax;

	//Subtract new label_edges from original label_edges
	imageCalculator("Subtract", label_edges, label_edges_edited);
	run("Red");
	setMinAndMax(0, 2);
	
	//Add nuclei outlines as overlay to the original image
	selectWindow(original);
	Overlay.remove;
	run("Add Image...", "image="+label_edges+" x=0 y=0 opacity=100 zero");
	run("Add Image...", "image="+label_edges_edited+" x=0 y=0 opacity=100 zero");
	close("label_edges");
//	draw_label_numbers(labelmap_edited, labelFontSize, fontColor);
	updateDisplay();
	close("label_edges");
	
	return newArray(labelmap_edited, label_edges_edited, nrNuclei);
}


function detectFoci(image, channel, fociSize, anisotropyFactor, firstTimeProcessing, labelmap_nuclei, labelmap_nuclei_3D, thresholdFactor) {

	//Determine threshold on the maximum projection (more robust), only where nuclei are present
	// Small FLAW: the area outside the nuclei is set to zero, but is included in the threshold calculation. (Could be fixed by setting 0 to NaN, but how to do that on the GPU?)
	// On the other hand: there is almost no difference, it seems.

	selectWindow(image);
	foci = "foci_ch"+channel;
	if(firstTimeProcessing == true) {
		rename(foci);
		Stack.setChannel(channel);
		Ext.CLIJ2_pushCurrentZStack(foci);
		rename(image);
		//Ext.CLIJ2_release(foci);
		if(debugMode) showImagefromGPU(foci);
	}
//	else if(debugMode) showImagefromGPU(foci);

	//Determine filter sizes
	if(fociSize == "tiny")	{ fociFilterSizeXY = 0;	fociSizeXY = 0.5;		}
	if(fociSize == "small")	{ fociFilterSizeXY = 0.5;	fociSizeXY = 1;		}
	if(fociSize == "average"){ fociFilterSizeXY = 1;	fociSizeXY = 1.5;	}
	if(fociSize == "large")	{ fociFilterSizeXY = 2;	fociSizeXY = 2;			}
	if(fociSize == "huge")	{ fociFilterSizeXY = 4;	fociSizeXY = 4;			}
	if(fociSize == "other") {
		fociSize = getNumber("Enter foci size in pixels", 2);
		fociFilterSizeXY = fociSize;
		fociSizeXY = fociSize;
	}

	Ext.CLIJ2_convertFloat(foci, foci_32bit);
	Ext.CLIJ2_pull(foci_32bit);
	rename(foci+"_outliers_removed");
	foci_outliers_removed = foci+"_outliers_removed";
	selectWindow(foci+"_outliers_removed");
	run("Remove Outliers", "block_radius_x="+minOf(25*fociSizeXY, 50)+" block_radius_y="+minOf(25*fociSizeXY, 50)+" standard_deviations=2 stack");
	if(debugMode) setBatchMode("show");
	Ext.CLIJ2_push(foci_outliers_removed);
//	if(gslices>1) Ext.CLIJ2_mask(foci_outliers_removed, labelmap_nuclei_3D, foci_outliers_removed_masked);
//	else Ext.CLIJ2_mask(foci_outliers_removed, labelmap_nuclei, foci_outliers_removed_masked);
//	showImagefromGPU(foci_outliers_removed_masked);
	//Determine the median of standard deviation signals in the oulier-removed nuclei
	run("Clear Results");
	Ext.CLIJ2_statisticsOfLabelledPixels(foci_outliers_removed, labelmap_nuclei);
	nucleiStddev_ = Table.getColumn("STANDARD_DEVIATION_INTENSITY", "Results");
	medianNucleiStddev = medianOfArray(nucleiStddev_);
	print("Median of the Standard Deviations of all nuclei: "+medianNucleiStddev);

	fociFilterSizeZ = 2*anisotropyFactor * fociFilterSizeXY;
	fociSizeZ = 2*anisotropyFactor * fociSizeXY;
	//Filter the foci - Difference of Gaussians, but then subtracting a blurred outlier-removed image
	if(gslices>1) Ext.CLIJ2_gaussianBlur3D(foci_32bit, foci_blurred, fociFilterSizeXY/2, fociFilterSizeXY/2, fociFilterSizeZ/2);
	else Ext.CLIJ2_gaussianBlur2D(foci_32bit, foci_blurred, fociFilterSizeXY/2, fociFilterSizeXY/2);
	if(gslices>1) Ext.CLIJ2_gaussianBlur3D(foci_outliers_removed, foci_outliers_removed_blurred, fociFilterSizeXY*2, fociFilterSizeXY*2, fociFilterSizeZ*2);
	else Ext.CLIJ2_gaussianBlur2D(foci_outliers_removed, foci_outliers_removed_blurred, fociFilterSizeXY*2, fociFilterSizeXY*2);

	if(debugMode) showImagefromGPU(foci_blurred);
	if(debugMode) showImagefromGPU(foci_outliers_removed_blurred);

	Ext.CLIJ2_release(foci_outliers_removed);
	foci_filtered = "foci_ch"+channel+"_filtered";
	Ext.CLIJ2_subtractImages(foci_blurred, foci_outliers_removed_blurred, foci_filtered);
	Ext.CLIJ2_release(foci_blurred);
	Ext.CLIJ2_release(foci_outliers_removed_blurred);

	if(debugMode) {
		print("Foci filter size XY: "+fociFilterSizeXY+", Z: "+fociFilterSizeZ);
		showImagefromGPU(foci_filtered);
		if(bits == 8) run("16-bit");
	}
//TO DO: For 3D, make max projection (see older versions)
//TO DO: Check Labelmap_detected_foci for 8-bit images. (no values >255)
//TO DO: Check fill holes (creates new labels in holes?)
	threshold = medianNucleiStddev * thresholdMultiplier;	//Set default threshold at n times the stddev (of the outlier-removed foci, which is a bit lower than the actual stddev, but less dependent on many foci being present)
	threshold = threshold*exp(thresholdFactor);
	print("Threshold default multiplier at "+thresholdMultiplier+", exponential bias set at "+thresholdFactor+": Threshold used: "+threshold);
	maxDisplaySetting = threshold * 10;

	if(foci_method == "AreaMaxima local maximum detection") {
		Ext.CLIJ2_pull(foci_filtered);
		selectWindow(foci_filtered);
//		setBatchMode("show");
//		hMin = 100;
		//setBatchMode(false);	//H-watershed does not work in batch mode if this image is not shown.
//		run("H_Watershed", "impin=["+foci+"] hmin="+hMin+" thresh="+threshold+" peakflooding="+H_max_peak_flooding + " outputmask=false allowsplitting=false");
//		print(getTitle);
		minSize = fociFilterSizeXY * fociFilterSizeXY * fociFilterSizeZ;
		run("AreaMaxima local maximum detection (2D, 3D)", "minimum="+minSize+" threshold="+threshold);
		rename("Labelmap_detected_foci_ch"+channel);
		Ext.CLIJ2_push("Labelmap_detected_foci_ch"+channel);
		close(foci_filtered);
		close("Labelmap_detected_foci_ch"+channel);
	}

	else if(foci_method == "Marker-controlled watershed") {
		//	Alternative function - could be faster: 
		//	Ext.CLIJx_detectAndLabelMaximaAboveThreshold(foci, detectedMaxima, 0, 0, 0, threshold, false)
		//	Ext.CLIJ2_pull(detectedMaxima);
		//	rename("detectedMaxima");

		//Detect all maxima
		allMaxima = "allMaxima";
		if(quasi2D && gslices>1) Ext.CLIJ2_detectMaximaSliceBySliceBox(foci_filtered, allMaxima, fociSizeXY, fociSizeXY);
		else Ext.CLIJ2_detectMaxima3DBox(foci_filtered, allMaxima, fociSizeXY, fociSizeXY, fociSizeZ);
		//remove spots lower than threshold
		MaskFociAboveThreshold = "MaskFociAboveThreshold_ch"+channel;
		Ext.CLIJ2_threshold(foci_filtered, MaskFociAboveThreshold, threshold);
		Ext.CLIJ2_mask(allMaxima, MaskFociAboveThreshold, maskedSpots);
		Ext.CLIJ2_release(allMaxima);
		if(debugMode) showImagefromGPU(MaskFociAboveThreshold);
		labeledSpots = "labeledSpots_ch"+channel;
		Ext.CLIJ2_connectedComponentsLabelingBox(maskedSpots, labeledSpots);	//Create label image with single pixels
		Ext.CLIJ2_release(maskedSpots);
//		if(bits == 8) Ext.CLIJ2_release(maskedSpots_16bit);
		if(debugMode) {
			showImagefromGPU(labeledSpots);
			//run("16-bit");
			run("glasbey on dark");
			run("Merge Channels...", "c1="+foci_filtered+" c2="+labeledSpots+" create keep");
			rename("Overlay_spots_ch"+channel);
			setBatchMode("show");
		}
//run("16-bit");
//run("glasbey on dark");
//run("Merge Channels...", "c1="+foci+" c2="+labeledSpots+" create keep");
//setBatchMode("show");
//Ext.CLIJx_detectAndLabelMaximaAboveThreshold(foci, detectedMaxima, 0, 0, 0, threshold, false)

		//Ext.CLIJx_seededWatershed(detectedMaxima, foci, Image_label_map_destination, threshold);	//Doesn't quite cut it - too many maxima detected?
		//Use CLIJx-MorphoLibJ marker-controlled watershed instead
		labelmap_foci_filled = "Labelmap_detected_foci_ch"+channel;
		//Invert foci before marker-controlled watershed!
		Ext.CLIJ2_invert(foci_filtered, foci_filtered_inverted);
		Ext.CLIJx_morphoLibJMarkerControlledWatershed(foci_filtered_inverted, labeledSpots, MaskFociAboveThreshold, labelmap_foci);
		Ext.CLIJx_morphoLibJFillHoles(labelmap_foci, labelmap_foci_filled);
		Ext.CLIJ2_release(labelmap_foci);

		Ext.CLIJ2_mask(labeledSpots, labelmap_nuclei_3D, labeledSpots_masked);
		Ext.CLIJ2_pullBinary(labeledSpots_masked);
		rename("foci_spots_ch"+channel);

		Ext.CLIJ2_release(foci_filtered);
		Ext.CLIJ2_release(foci_filtered_inverted);
		Ext.CLIJ2_release(labeledSpots);
		Ext.CLIJ2_release(labeledSpots_masked);
		Ext.CLIJ2_release(MaskFociAboveThreshold);
	}
	labelmapDetectedFoci = "Labelmap_detected_foci_ch"+channel;

	//Exclude foci too small or too large
	if(debugMode) {
		showImagefromGPU(labelmapDetectedFoci);
		run("glasbey on dark");
		resetMinAndMax;
	}
	Ext.CLIJ2_excludeLabelsOutsideSizeRange(labelmapDetectedFoci, labelmapDetectedFociFiltered, minFociSize, maxFociSize);
	if(debugMode) showImagefromGPU(labelmapDetectedFoci);
	//Exclude foci outside of nuclei
	Ext.CLIJ2_mask(labelmapDetectedFociFiltered, labelmap_nuclei_3D, labelmapDetectedFociFilteredInsideNuclei);
	Ext.CLIJ2_release(labelmapDetectedFociFiltered);
	Ext.CLIJ2_closeIndexGapsInLabelMap(labelmapDetectedFociFilteredInsideNuclei, labelmap_foci_final);
	Ext.CLIJ2_release(labelmapDetectedFociFilteredInsideNuclei);
	Ext.CLIJ2_getMaximumOfAllPixels(labelmap_foci_final, max);
	print("\n"+max+" foci detected in "+nrNuclei+" nuclei in channel "+channel + "\n");
	Ext.CLIJ2_pull(labelmap_foci_final);
	if(debugMode) {
		run("glasbey on dark");
		resetMinAndMax;
		setBatchMode("show");
	}
	rename("Labelmap_detected_foci_filtered_ch"+channel);
//setMinAndMax(0, max);
//run("glasbey_on_dark");
//setBatchMode("show");

	//create binary from the labelmap
	Ext.CLIJ2_pullBinary(labelmap_foci_final);
	rename("Mask_foci_ch"+channel);

	Ext.CLIJ2_release(labelmap_foci_final);
	if(!debugMode) close(foci+"_outliers_removed");

	return "Labelmap_detected_foci_filtered_ch"+channel;
}


function createLabelmapNuclei3D(labelmap_nuclei) {	//Create 3D labelmap image by stacking the 2D labelmap 
//	Ext.CLIJ2_push(labelmap_nuclei);
//	labelmap_nuclei = "Labelmap_nuclei";
	labelmap_nuclei_3D = "Labelmap_nuclei_3D";
	selectWindow("nuclei");
	labelmap_nuclei_3D = "nuclei";
	Ext.CLIJ2_push(labelmap_nuclei_3D);	//Push this image to the GPU and turn it into a labelmap by copying slices. Another way is to create it from scratch in the GPU and then fill it.
	if(gslices > 1) {
		for (i = 0; i < gslices; i++) {
			Ext.CLIJ2_copySlice(labelmap_nuclei, labelmap_nuclei_3D, i);
		}
	}
	else labelmap_nuclei_3D = labelmap_nuclei;	//For compatibility with 2D images
	return labelmap_nuclei_3D;
}


function mergeOriginalAndDetection(original, nrNuclei, nuclei_outlines, foci_mask, foci_spots, channel) {
	// Merge the original foci channel with the detected foci
	selectWindow("nuclei");
	run(bits+"-bit");	//convert to allow merging
	foci_RAW = "foci_ch"+channel;
	if(bits == 16) Ext.CLIJ2_convertUInt16(foci_RAW, foci_RAW_8or16bit);
	else if(bits == 8) Ext.CLIJ2_convertUInt8(foci_RAW, foci_RAW_8or16bit);
	Ext.CLIJ2_pull(foci_RAW_8or16bit);	//Still in GPU memory
	if(optimizationMode == true && doneOptimizing == true) Ext.CLIJ2_release(foci_RAW);
	Ext.CLIJ2_release(foci_RAW_8or16bit);

	selectWindow(foci_mask);
	run(bits+"-bit");	//convert to allow merging
	selectWindow(foci_spots);
	run(bits+"-bit");	//convert to allow merging
	run("Merge Channels...", "c1=nuclei c2=" + foci_RAW_8or16bit + " c3=" + foci_mask + " c4=" + foci_spots + " create keep");
	getDimensions(width, height, channels, slices, frames);
	Stack.setSlice(slices/2);
	Stack.setChannel(1);
    setLut(reds, greens, blues);
 	run("Enhance Contrast...", "saturated=0.02");  
	Stack.setChannel(2);
	run("Green");
	run("Enhance Contrast...", "saturated=0.02");
	getMinAndMax(min, max);
	setMinAndMax(min, maxOf(max, min + maxDisplaySetting));
	Stack.setChannel(3);
	resetMinAndMax;
	run("Magenta");
	Stack.setChannel(4);
	resetMinAndMax;
	run("Grays");	
	rename("Foci_overlay_ch"+channel);

	//Add nuclei outlines as overlay to the original image
	selectWindow("Foci_overlay_ch"+channel);
	for (i = 1; i <= gslices; i++) {
		if(gslices>1) Stack.setSlice(i);
		run("Add Image...", "image="+nuclei_outlines+" x=0 y=0 opacity="+labelOpacity+" zero");	//Add labelmap to image as overlay
		Overlay.setPosition(0, 0, 0);
	}
	if(addNumbersOverlay) {
		run("Clear Results");
		Ext.CLIJ2_statisticsOfLabelledPixels(labelmap_nuclei, labelmap_nuclei); //Somehow if you put (image, labelmap) as arguments the pixel count is wrong
		setFont("SansSerif", labelFontSize, "bold antialiased");
		color = color_to_hex(fontColor);
		setColor(color);
		for (i = 0; i < nrNuclei; i++) {
			x = getResult("MASS_CENTER_X", i);
			y = getResult("MASS_CENTER_Y", i);
			Overlay.drawString(i+1, x - labelFontSize/2, y + labelFontSize/2);
		}
	}
}



function measureFoci(original, channel, nrNuclei, labelmap_nuclei_3D, labelmap_foci, boundingBox_X, boundingBox_Y, boundingBox_width, boundingBox_height) {
	Ext.CLIJ2_push(labelmap_foci);
	if(debugMode) showImagefromGPU(labelmap_foci);
	selectWindow(original);
	run("Duplicate...", "title=foci_RAW_ch"+channel + " duplicate channels=" + channel);

	nucleus_id_ = newArray(nrNuclei);
	nucleus_area_ =  newArray(nrNuclei);
	foci_count_ = newArray(nrNuclei);
	foci_mean_int_ = newArray(nrNuclei);
	foci_median_int_ = newArray(nrNuclei);
	foci_size_ = newArray(nrNuclei);
	foci_median_size_ = newArray(nrNuclei);
	foci_sum_ = newArray(nrNuclei);
	nucleus_mean_int_ = newArray(nrNuclei);
	background_ = newArray(nrNuclei);

	//Loop over all nuclei
	startTime = getTime();
	foci_raw = "foci_RAW_ch"+channel;
	Ext.CLIJ2_push(foci_raw);
	for (i = 0; i < nrNuclei; i++) {
		showStatus("Analyzing foci in channel "+channel+", nucleus "+i+1 + "/" + nrNuclei);
		showProgress(i, nrNuclei);
		nrResults = nResults;
		labelmap_foci_nucleus = "labelmap_foci_nucleus";

		// statisticsOfLabelledPixels and closeIndexGapsInLabelMap are very slow for large images.
		// Crop labelmaps and raw images before measuring - much faster for large images
		Ext.CLIJ2_crop3D(foci_raw, foci_raw_cropped, boundingBox_X[i], boundingBox_Y[i], 0, boundingBox_width, boundingBox_height, gslices);
		if(debugMode) {
//			showImagefromGPU(foci_raw_cropped);
//			getDimensions(width, height, channels, slices, frames);
//			print("dimensions of nucleus "+i+1+" image: "+width+ ", "+height);
		}
		Ext.CLIJ2_crop3D(labelmap_nuclei_3D, labelmap_nuclei_cropped, boundingBox_X[i], boundingBox_Y[i], 0, boundingBox_width, boundingBox_height, gslices);
		Ext.CLIJ2_crop3D(labelmap_foci, labelmap_foci_cropped, boundingBox_X[i], boundingBox_Y[i], 0, boundingBox_width, boundingBox_height, gslices);
		Ext.CLIJ2_maskLabel(labelmap_foci_cropped, labelmap_nuclei_cropped, labelmap_foci_nucleus, i+1);

		//This is roughly equally fast as maskLabel:
		//Ext.CLIJ2_labelToMask(labelmap_nuclei_3D, mask_nucleus, i+1);
		//Ext.CLIJ2_mask(labelmap_foci, mask_nucleus, labelmap_foci_nucleus);
		//Also try:
		//Ext.CLIJ2_maskStackWithPlane(Image_source, Image_mask, Image_destination);

		Ext.CLIJ2_closeIndexGapsInLabelMap(labelmap_foci_nucleus, labelmap_foci_nucleus_relabeled);
	
		Ext.CLIJ2_statisticsOfLabelledPixels(foci_raw_cropped, labelmap_foci_nucleus_relabeled);
		
		//Count the number of foci
		Ext.CLIJ2_getMaximumOfAllPixels(labelmap_foci_nucleus_relabeled, nrFoci);
		foci_count_[i] = nrFoci;

		//Get relevant results and store in arrays
		mean_ = newArray(nrFoci);
		size_ = newArray(nrFoci);
		sum = 0;
		for(k = 0; k < nrFoci; k++) {
			setResult("Nucleus", nrResults + k, i+1);
			mean_[k] = getResult("MEAN_INTENSITY", nrResults + k);
			if(gslices > 1) size_[k] = getResult("PIXEL_COUNT", nrResults + k) * pixelWidth * pixelHeight * pixelDepth;
			else size_[k] = getResult("PIXEL_COUNT", nrResults + k) * pixelWidth * pixelHeight;
			sum = sum + getResult("SUM_INTENSITY", nrResults + k);
		}
		Array.getStatistics(mean_, min, max, meanIntensity, stdDev);
		if(mean_.length > 0) medianIntensity = medianOfArray(mean_);
		else medianIntensity = NaN;
		Array.getStatistics(size_, min, max, meanSize, stdDev);
		if(mean_.length > 0) medianSize = medianOfArray(size_);
		else medianSize = NaN;
		
		nucleus_id_[i] = i+1;
		foci_count_[i] = nrFoci;
		foci_mean_int_[i] = meanIntensity;
		foci_median_int_[i] = medianIntensity;
		foci_size_[i] = meanSize;
		foci_median_size_[i] = medianSize;
		foci_sum_[i] = sum;
		
		if(i%100 == 0) run("Clear Results");	//A large Results table is slowing down a lot! Intermediate clearing helps.
	}
	Ext.CLIJ2_release(foci_raw_cropped);
	Ext.CLIJ2_release(labelmap_nuclei_cropped);
	Ext.CLIJ2_release(labelmap_foci_cropped);
	Ext.CLIJ2_release(labelmap_foci_nucleus_relabeled);
	Ext.CLIJ2_release(labelmap_foci_nucleus);
	print("Measuring all foci in channel "+channel+" took "+ getTime() - startTime +" ms.");

	//Measure overall nuclear intensity
	run("Clear Results");
	Ext.CLIJ2_statisticsOfBackgroundAndLabelledPixels(foci_raw, labelmap_nuclei_3D);
	background = getResult("MEAN_INTENSITY", 0);
	Array.fill(background_, background);
	for (i = 0; i < nResults-1; i++) {
		nucleus_mean_int_[i] = getResult("MEAN_INTENSITY", i+1);
	}

	Ext.CLIJ2_release(foci_raw);

	//Add relevant results to the resultTable
	Table.setColumn("Background intensity ch"+channel, background_, resultTable);
	Table.setColumn("Mean nucleus intensity ch"+channel, nucleus_mean_int_, resultTable);
	Table.setColumn("Foci count ch"+channel, foci_count_, resultTable);
	Table.setColumn("Mean foci intensity ch"+channel, foci_mean_int_, resultTable);
	Table.setColumn("Median foci intensity ch"+channel, foci_median_int_, resultTable);
	if(gslices > 1) Table.setColumn("Mean foci volume ("+unit+"^3) ch"+channel, foci_size_, resultTable);
	else Table.setColumn("Mean foci area ("+unit+"^2) ch"+channel, foci_size_, resultTable);
	if(gslices > 1) Table.setColumn("Median foci volume ("+unit+"^3) ch"+channel, foci_median_size_, resultTable);
	else Table.setColumn("Median foci area ("+unit+"^2) ch"+channel, foci_median_size_, resultTable);
	Table.setColumn("Total foci intensity ch"+channel, foci_sum_, resultTable);
	Table.update();
}


function measure_nuclear_intensities(original, nrNuclei, labelmap_nuclei_3D, gchannels, fociChannelA, fociChannelB, resultTable) {
	background_ = newArray(nrNuclei);
	for(c=1; c<=gchannels; c++) {
		if(c != fociChannelA && c != fociChannelB) {
			selectWindow(original);
			run("Clear Results");
			Stack.setChannel(c);
			channel_to_measure = original+"_ch"+c;
			run("Duplicate...", "title=[" + channel_to_measure + "] duplicate channels=" + c);
			Ext.CLIJ2_push(channel_to_measure);
			Ext.CLIJ2_statisticsOfBackgroundAndLabelledPixels(channel_to_measure, labelmap_nuclei_3D);
			intensity_ = Table.getColumn("MEAN_INTENSITY", "Results");
			background = intensity_[0];
			Array.fill(background_, background);
			intensity_ = Array.slice(intensity_, 1, intensity_.length);
			Table.setColumn("Background intensity ch"+c, background_, resultTable);
			Table.setColumn("Mean nucleus intensity ch"+c, intensity_, resultTable);
			Table.update;
			Ext.CLIJ2_release(channel_to_measure);
		}
	}
}


function computeOverlap(mask_fociA, mask_fociB, nrNuclei, labelmap_nuclei_3D, boundingBox_X, boundingBox_Y, boundingBox_width, boundingBox_height) {
	showStatus("Computing foci colocalization...");
	startTime = getTime();
	run("Clear Results");
	Ext.CLIJ2_push(mask_fociA);
	Ext.CLIJ2_push(mask_fociB);
	Ext.CLIJ2_binaryIntersection(mask_fociA, mask_fociB, overlap);
	Ext.CLIJ2_connectedComponentsLabelingBox(overlap, labelmap_overlap);
//Ext.CLIJ2_generateBinaryOverlapMatrix(labelmap_fociA, labelmap_fociB, overlap_matrix);
//Ext.CLIJ2_pull(overlap_matrix);
//setBatchMode("show");
	Ext.CLIJ2_getJaccardIndex(mask_fociA, mask_fociB, jaccardIndex);
	//print("Jaccard Index (before overlap size exclusion): "+jaccardIndex);

	Ext.CLIJ2_release(mask_fociA);
	Ext.CLIJ2_release(mask_fociB);
	Ext.CLIJ2_release(overlap);

	//Filter on overlap size
	labelmap_overlap_filtered = "labelmap_overlap_filtered";
	Ext.CLIJ2_excludeLabelsOutsideSizeRange(labelmap_overlap, labelmap_overlap_filtered, minOverlapSize, maxFociSize);
	Ext.CLIJ2_release(labelmap_overlap);
	Ext.CLIJ2_getMaximumOfAllPixels(labelmap_overlap_filtered, max);
	print(max+" colocalizing foci detected with overlap > "+ minOverlapSize +" voxels.");

	//Count overlapping foci for all nuclei
	overlap_count_ = newArray(nrNuclei);
	for (i = 0; i < nrNuclei; i++) {
		showStatus("Computing foci overlap in nucleus "+i+1 + "/" + nrNuclei);
		showProgress(i, nrNuclei);

		labelmap_overlap_nucleus = "labelmap_overlap_nucleus_"+i+1;

		// Crop labelmaps (nuclei and overlap) before measuring - much faster for large images
		Ext.CLIJ2_crop3D(labelmap_nuclei_3D, labelmap_nuclei_cropped, boundingBox_X[i], boundingBox_Y[i], 0, boundingBox_width, boundingBox_height, gslices);
		Ext.CLIJ2_crop3D(labelmap_overlap_filtered, labelmap_overlap_cropped, boundingBox_X[i], boundingBox_Y[i], 0, boundingBox_width, boundingBox_height, gslices);
		if(debugMode) {
//			showImagefromGPU(labelmap_nuclei_cropped);
//			getDimensions(width, height, channels, slices, frames);
//			print("dimensions of nucleus "+i+1+" image: "+width+ ", "+height);
		}
		Ext.CLIJ2_maskLabel(labelmap_overlap_cropped, labelmap_nuclei_cropped, labelmap_overlap_nucleus, i+1);

		Ext.CLIJ2_closeIndexGapsInLabelMap(labelmap_overlap_nucleus, labelmap_overlap_nucleus_relabeled);
		Ext.CLIJ2_getMaximumOfAllPixels(labelmap_overlap_nucleus_relabeled, nrOverlappingFoci);
		
		overlap_count_[i] = nrOverlappingFoci;

		//Cannot re-use the memory, because the shape is different for every nucleus!
		Ext.CLIJ2_release(labelmap_nuclei_cropped);
		Ext.CLIJ2_release(labelmap_overlap_cropped);
		Ext.CLIJ2_release(labelmap_overlap_nucleus);
		Ext.CLIJ2_release(labelmap_overlap_nucleus_relabeled);
	}
	Ext.CLIJ2_release(labelmap_overlap_filtered);

	//Add overlap counts to the resultTable
	Table.setColumn("Overlapping foci", overlap_count_, resultTable);
	fociCountChannelA_ = Table.getColumn("Foci count ch"+fociChannelA, resultTable);
	overlap_fraction_ = divideArrays(overlap_count_, fociCountChannelA_);
	overlap_percentage_ = multiplyArraywithScalar(overlap_fraction_, 100);
	Table.setColumn("Overlap percentage (ch"+fociChannelA+" -> ch"+fociChannelB+")", overlap_percentage_, resultTable);

// TO DO: Use CLIJx_labelOverlapCountMap to calculate the number of overlapping foci.
// It will replace the above code, except for the overlap size filtering, which may not be possible.
// Also, you can use it to count foci in every cell by passing the nuclei labelmap and a foci labelmap.

	//Create Results table column for the overlap visualization, including the background (entry 0)
	run("Clear Results");
	setResult("Overlap percentage", 0, 0);
	for (i = 0; i < nrNuclei; i++) {
		setResult("Overlap percentage", i+1, overlap_percentage_[i]);
	}
//	updateResults();

	//Visualize the amount of overlapping foci in a map
	//The following line doesn't work. Why not?
//	Ext.CLIJ2_generateParametricImageFromResultsTableColumn(labelmap_nuclei_3D, overlap_count_map, "Overlap percentage");

	Ext.CLIJ2_pushResultsTableColumn(overlap_fill, "Overlap percentage");
		if(debugMode == true) showImagefromGPU(overlap_fill);
	overlap_count_map = "overlap_count_map";
	Ext.CLIJ2_replaceIntensities(labelmap_nuclei_3D, overlap_fill, overlap_count_map);
	Ext.CLIJ2_release(overlap_fill);
		if(debugMode == true) showImagefromGPU(overlap_count_map);
	if(bits == 16) Ext.CLIJ2_convertUInt16(overlap_count_map, overlap_count_map_8or16bit);
	else if(bits == 8) Ext.CLIJ2_convertUInt8(overlap_count_map, overlap_count_map_8or16bit);
	Ext.CLIJ2_release(overlap_count_map);
	overlap_fraction_map = "overlap_fraction_map";
	if(bits == 16) Ext.CLIJ2_replaceIntensity(overlap_count_map_8or16bit, overlap_fraction_map, 65535, 0);
	else if(bits == 8) Ext.CLIJ2_replaceIntensity(overlap_count_map_8or16bit, overlap_fraction_map, 255, 0);
	Ext.CLIJ2_pull(overlap_fraction_map);
	Ext.CLIJ2_release(overlap_fraction_map);
	Ext.CLIJ2_release(overlap_count_map_8or16bit);

	print("Measuring colocalization took "+ getTime() - startTime +" ms.");

	//Create foci mask overlap image 
	run("Merge Channels...", "c1="+mask_fociA+" c2="+mask_fociB+" c3="+overlap_fraction_map+" create");
	Stack.setChannel(1);
	run("Green");
	Stack.setChannel(2);
	run("Red");
	rename("Foci overlap map");
	Stack.setChannel(3);
	setLut(reds, greens, blues);
	setMinAndMax(0, 100);
	setBatchMode("show");

	//Add nuclei outlines as overlay to the overlap image
	selectWindow("Foci overlap map");
	for (i = 1; i <= gslices; i++) {
		if(gslices > 1) Stack.setSlice(i);
		run("Add Image...", "image="+nuclei_outlines+" x=0 y=0 opacity="+labelOpacity+" zero");	//Add labelmap to image as overlay
	}
	if(addNumbersOverlay) {
		run("Clear Results");
		Ext.CLIJ2_statisticsOfLabelledPixels(labelmap_nuclei, labelmap_nuclei); //Somehow if you put (image, labelmap) as arguments the pixel count is wrong
		setFont("SansSerif", labelFontSize, "bold antialiased");
		color = color_to_hex(fontColor);
		setColor(color);
		for (i = 0; i < nrNuclei; i++) {
			x = getResult("MASS_CENTER_X", i);
			y = getResult("MASS_CENTER_Y", i);
			Overlay.drawString(i+1, x - labelFontSize/2, y + labelFontSize/2);
		}
	}
	if(gslices > 1) Stack.setSlice(gslices/2);
	setBatchMode("show");

}

//Create a nice azure blue LUT for the nuclei overlay
function create_azure_lut(reds, greens, blues) {
	for (i=0; i<256; i++) {
		reds[i] = 0;
		greens[i] = 0.5*i;
		blues[i] = i;
	}
}

//Convert a color into a hexadecimal code
function color_to_hex(color) {
	colorArray = split(color,",,");
	hexcolor = "#" + IJ.pad(toHex(colorArray[0]),2) + IJ.pad(toHex(colorArray[1]),2) + IJ.pad(toHex(colorArray[2]),2);
	return hexcolor;
}


//Returns the median of an array
function medianOfArray(array) {
	Array.sort(array);
	return(array[floor(array.length/2)]);
}

//Divides the elements of two arrays and returns the new array
function divideArrays(array1, array2) {
	divArray=newArray(lengthOf(array1));
	for (a=0; a<lengthOf(array1); a++) {
		divArray[a]=array1[a]/array2[a];
	}
	return divArray;
}

//Multiplies all elements of an array with a scalar and returns the new array
function multiplyArraywithScalar(array, scalar) {
	multiplied_array=newArray(lengthOf(array));
	for (a=0; a<lengthOf(array); a++) {
		multiplied_array[a]=array[a]*scalar;
	}
	return multiplied_array;
}

function showImagefromGPU(string_Image) {
	Ext.CLIJ2_pull(string_Image);
	setBatchMode("show");
}

function cleanup() {
//	close("nuclei");
	close("foci");
}

