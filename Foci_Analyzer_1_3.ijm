#@ File[] (label = "Input files", style="File") files
#@ File (label = "Output folder", style = "directory") outputFolder

#@ Integer (label = "Nuclei channel", value = 1) nucleiChannel
#@ Integer (label = "Foci channel A", value = 2) fociChannelA
#@ Integer (label = "Foci channel B", value = 3) fociChannelB
#@ Boolean (label = "Also detect foci in channel B and perform colocalization?", persistence=true, value=true, description="Check this box to calculate colocalization between two foci channels.") detectFociChannelB

#@ String (label = "3D image handling", choices={"Analyze foci in 3D", "Detect foci on the Maximum Intensity Projection", "Use quasi-2D foci detection (detect foci separately in every Z-slice)", "Process a single z-slice only (specify below which slice)"}, style="radioButtonVertical") ThreeDHandling
#@ Integer (label = "[single z-slice foci detection only] Slice nr", value=1, min=1) singleSlice
#@ Integer (label = "Remove image borders (pixels)", value = 0, min=0) cropBorder
#@ Integer (label = "Image XY binning before analysis [1-n]", value = 1, min=1, description="Image binning will greatly shorten the analysis time") XYBinning

#@ String (value="Nuclei detection settings", visibility="MESSAGE") nuclei_message
#@ String (label = "Nuclei segmentation method", choices={"StarDist nuclei segmentation (recommended)","Cellpose cytoplasm segmentation", "Classic nuclei segmentation"}, style="listBox") nucleiSegmentationChoice
#@ Integer (label = "Stardist nuclei binning factor [1-n]", value = 1, min=1, description="A binning >1 for high-resolution input images can prevent 'oversegmentation', and will also somewhat speed up processing.") downsampleFactorStarDist
#@ Double (label = "Probablility threshold [0.0-1.0] (StarDist/Cellpose)", value = 0.5, min=0, max=1, description="lower values will accept more nuclei") probabilityThreshold

#@ Integer (label = "Remove nulei with diameter smaller than (units)", value = 4) minNucleusSize_setting
#@ Integer (label = "Remove nulei with diameter larger than (units)", value = 50) maxNucleusSize_setting
#@ Boolean (label = "Exclude nulei on image edges", value = true) excludeOnEdges
#@ String (label = "Manual nuclei removal", choices={"No thanks","Manually remove nuclei", "load previously saved removals (from output folder)"}, value = "No thanks", description="Allow nuclei removal by clicking, or load previous removals (Masks should be in the output folder)") manualRemoveNuclei

#@ String (value="Foci detection settings", visibility="MESSAGE") foci_message1
#@ Boolean (label = "Enable foci parameters optimization mode?", value=true) optimizationModeSetting
#@ String (label = "Foci size channel A", choices={"tiny","small","average","large","huge","other"}, style="listBox", value="average", description="Simplified foci detection parameter") fociSizeA
#@ String (label = "Foci size channel B", choices={"tiny","small","average","large","huge","other"}, style="listBox", value="average", description="Simplified foci detection parameter") fociSizeB
#@ String (label = "Foci detection method", choices={"Marker-controlled watershed (recommended)","AreaMaxima local maximum detection"}, style="listBox") foci_method

#@ Double (label = "Foci intensity threshold bias channel A", value = 0.0, min=-2.5, max=2.5, style="scroll bar", stepSize=0.01) thresholdFactorA
#@ Double (label = "Foci intensity threshold bias channel B", value = 0.0, min=-2.5, max=2.5, style="scroll bar", stepSize=0.01) thresholdFactorB
#@ Double (label = "Minimum foci size (pixels/voxels)", value = 3) minFociSize
#@ Double (label = "Maximum foci size (pixels/voxels)", value = 9999) maxFociSize

#@ Integer (label = "Minimum overlap of foci to colocalize (pixels/voxels)", min = 1, value = 1, description="Foci in channel A and B will be counted as colocalizing only if they overlap with at least this area/volume") minOverlapSize

#@ String (value="Visualization options", visibility="MESSAGE") coloc_message
#@ String(label = "Nuclei outline color", choices={"White","Red","Green","Blue","Cyan","Magenta","Yellow"}, value = "Cyan") outlineColor
#@ ColorRGB(label = "Nuclei label color", value="orange") fontColor
#@ Integer (label = "Nuclei label size", value=12, min=4) labelFontSize

#@ Boolean (label = "Debug mode (show intermediate images)", value=false) debugMode

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
 * Author: Bram van den Broek, The Netherlands Cancer Institute, b.vd.broek@nki.nl
 * 
 * 
 * 
 * Changelog 
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
 * Version 1.0, July 2022:
 * - Added possibility to manually remove segmented nuclei by clicking. Masks are saved and can be loaded for re-analysis.
 * 
 * Version 1.1, July 2022:
 * - Fixed a small bug: nuclear intensity in foci channel B is now measured also when detectFociChannelB is false.
 * - Changes in AreaMaximum detection (minimum foci size) 
 * - Now asks to create output folder if it doesn't exist.
 * - Classic nuclei segmentation works again
 * - Fixed a bug where the nucleus area was not correctly reported when Stardist nuclei binning factor > 1
 * 
 * Version 1.2, October 2022:
 * - Fixed critical bug where after threshold optimization of the second image Fiji would crash (because an image was released from the GPU)
 * - Improved log window readability
 * - Fixed a bug (in debug mode) where Fiji crashed after threshold optimization if the overlay image was not selected
 * 
 * Version 1.22, October 2022:
 * - Fixed a bug where manual nuclei removal didn't work on max projections
 * - Added labels to channels in the overlay image
 * - Implemented a maximum tile size for Stardist
 * 
 * Version 1.23, October 2022:
 * - Added option to crop image borders
 * 
 * Version 1.24, December 2022:
 * - Added total nucleus intensity (mean * area) as output column
 * - Save cropBorders and XYBinning in the output image metadata
 * - check if Foci channel B actually exists (before it would just take the last channel) and is not the same as Foci channel A
 * 
 * Version 1.25, December 2022
 * - Fixed bug where foci spots were not displayed correctly when threshold was changed
 * - Improved visualization options
 * 
 * Version 1.3, December 2022
 * - If 'Alt' is pressed before opening a file, the user can select a ROI to do a preview analysis on.
 * - downscaled ROIs are now smoothed using spline fitting
 * - Improved GUI: 3D image handling radiobuttons instead of separate checkboxes
 * - Improved logging: Threshold factor * bias is now mentioned in the log window
 * - Colocalization overlap image is now also 3D and saved as .Gif file 
 * 
*/

version = 1.3;

/* TODO / IDEAS
 *  
 * Fix channel 4 (foci spots) for filtered foci. Currently, max foci size does not remove the spots.
 * Do not remove lines in log window when areaMaximum detection is used.
 *  
 * Use ALT to pop up options menu (background subtraction, nuclear segmentation, foci size, etc.
 * and allow the user to determine the right settings! See MosaicExplorerJ macro.
 * Currently, pressing ALT during opening the image allows for making a selection.
 * 
 * Upscale and smooth ROIs (spline interpolate met Roi.setPolygonSplineAnchors(x, y)?)
 * 
 * Add more: description="..." to the options (onmouseover)
 * 
 * Open and split timelapse image; save data as individual .tsv files? (in order to prevent too much bookkeeping)
 * 
 * Background intensity is now very crude (just the mean outside the nuclei). So: mask the background, set (percentile) threshold, measure mean/median. 
 * 
 * Auto-bin nuclei based on pixel size
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

maxTileSize = 2048;			// Maximum StarDist tile size
thresholdMultiplier = 3;	//Default threshold multiplier - foci should be this number of times more intense than the standard deviation of the background (average of all nuclei) 
minFociSize = minFociSize;
maxFociSize = maxFociSize;

H_max_peak_flooding = 100;	//Percentage of the peak
if(outlineColor == "White") outlineColor = "Grays";

thickOutlines = false;		//Width of nuclei outlines, 1 or 2 pixels
thickOutlinesThresholdSize = 1200;	//Above this image size nuclei outlines are thick

labelOpacity = 100;			//opacity of outlines
addNumbersOverlay = true;
saveOverlayImage = true;

var optimizationMode = optimizationModeSetting;
var processAllOtherImages = false;
var useROI = false;
var reAnalyzeFullImage = false;
var currentFileNr = 0;
var doneOptimizing;
var displayMode = "composite";
var activeChannels = "1111";
var maxDisplaySetting;
var overlayBrightness = "bright";
var processTime = 0;
var threshold = 0;


//Create the azure LUT for the nuclei (TO DO: do this in a more elegant way than globals)
reds = newArray(256);
greens = newArray(256); 
blues = newArray(256);
create_azure_lut(reds, greens, blues);
//create_lut(reds, greens, blues)

saveSettings();

run("Set Measurements...", "area mean standard integrated median redirect=None decimal=3");
run("Conversions...", " ");
setOption("BlackBackground", true);
setForegroundColor(255, 255, 255);
run("Colors...", "foreground=white background=black selection=cyan");

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
if(!File.exists(outputFolder)) {
	createOutputFolder = getBoolean("The output folder '"+outputFolder+"' does not exist. Create?", "Of course, go ahead!", "See if I care!");
	if(createOutputFolder) File.makeDirectory(outputFolder);
	else {
		formatHardDrive = getBoolean("Think you're funny, eh? Try this:\nFormat the hard drive?", "Yes, goodbye forever", "No! Please mr. Foci Analyzer, I'll do anything you ask!");
		if(formatHardDrive) {
			showMessage("Ok, you wished for it!");
			exit("Oh, wait. I'll erase myself as well.\nCall it your lucky day then!");
		}
		else exit("That's more like it. Now, run the macro and try again.");
		if(createOutputFolder) File.makeDirectory(outputFolder);
	}
}
for (currentFileNr = 0; currentFileNr < nrOfImages; currentFileNr++) {
	print("\nProcessing file "+currentFileNr+1+"/"+nrOfImages+": "+files[currentFileNr] + "\n");
	processFile(currentFileNr, files[currentFileNr], outputFolder);
}
print("\n-------------------------------------------------------------------------");
print("Finished processing "+nrOfImages+" files in "+processTime*60+" seconds ("+d2s(processTime,1)+" minutes).");
print("Average speed: "+d2s((nrOfImages)/processTime,1)+" images per minute.");
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
		if(reAnalyzeFullImage == false && optimizationModeSetting == true && processAllOtherImages == false) optimizationMode = true;
		process_current_series(file, true);
		if(useROI) {
			reAnalyzeFullImage = true;
			currentFileNr--;	//ugly change of global variable, but I can't see another easy way
			if(optimizationModeSetting == true && processAllOtherImages == false) optimizationMode = false;
		}
		else reAnalyzeFullImage = false;
	}
	else {	//Use Bio-Formats
		for(s = 0; s < nr_series; s++) {
			run("Close All");
			run("Bio-Formats Importer", "open=["+file+"] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT series_"+s+1);
			seriesName = getTitle();
			seriesName = replace(seriesName,"\\/","-");	//replace slashes by dashes in the seriesName
	//		print(seriesName);
	//		outputPath = output + File.separator + substring(seriesName)

			if(XYBinning > 1) run("Bin...", "x="+XYBinning+" y="+XYBinning+" z=1 bin=Average");
			process_current_series(seriesName, false);
			if(useROI) s--;	//Analyze this series again without ROI
			//may work, but only for files opened by bioformats
		}
	}
}



function process_current_series(image, nameHasExtension) {
	if(ThreeDHandling == "Process a single z-slice only (specify below which slice)") run("Duplicate...", "duplicate slices="+singleSlice);
	Stack.setDisplayMode("grayscale");
	Stack.setChannel(nucleiChannel);
	original = getTitle();

	if(nameHasExtension) imageName = File.getNameWithoutExtension(image);
	else imageName = image;
	
	//Initialize image and table
	getDimensions(gwidth, gheight, gchannels, gslices, gframes); // global variables
	getVoxelSize(pixelWidth, pixelHeight, pixelDepth, unit);
	if(cropBorder>0) {
		makeRectangle(cropBorder, cropBorder, gwidth-2*cropBorder, gheight-2*cropBorder);
		run("Crop");
		getDimensions(gwidth, gheight, gchannels, gslices, gframes); // global variables
	}
	if(gwidth > thickOutlinesThresholdSize && gheight > thickOutlinesThresholdSize) {
		thickOutlines = true;
	}
	
	if(gslices > 1) anisotropyFactor = pixelWidth / pixelDepth;
	else anisotropyFactor = 0;
	bits = bitDepth();
	
//TO DO: FIND A BETTER WAY FOR THIS:
	if(bits == 8) run("16-bit");	//Convert to 16-bit, because 8-bit restricts the foci labelmap to 255 foci

	Stack.setSlice(gslices/2);
	setBatchMode("show");
	run("Enhance Contrast", "saturated=0.35");

	resultTable = "Foci results";
	if(isOpen(resultTable)) close(resultTable);
	Table.create(resultTable);
	Table.setLocationAndSize(0, 0, 1000, 500);

	if(gslices > 1 && ThreeDHandling == "Detect foci on the Maximum Intensity Projection") {
		selectWindow(original);
		run("Z Project...", "projection=[Max Intensity]");
		setBatchMode("show");
		original = getTitle();
		gslices = 1;
	}

	//Run analysis only on a ROI
	wait(50);
	useROI = false;
	if(isKeyDown("alt")) {
		setTool("rectangle");
		waitForUser("Select a ROI to analyze");
		if(selectionType==0) {
			useROI = true;
			getSelectionBounds(x_ROI, y_ROI, width_ROI, height_ROI);
			run("Duplicate...", "duplicate title=selection");
			process_image = getTitle();
		}
		else {
			print("Ignoring non-rectangular selection...");
			run("Select None");
		}
	}
	if(useROI == false) {
		process_image = original;
		x_ROI = 0;
		y_ROI = 0;
	}

	//Segment and label nuclei
	if (nucleiSegmentationChoice == "StarDist nuclei segmentation (recommended)") nuclei_info = segmentNucleiStarDist(process_image, nucleiChannel, probabilityThreshold, pixelWidth, unit, resultTable);
	else if(nucleiSegmentationChoice == "Classic nuclei segmentation") nuclei_info = segmentNucleiClassic(process_image, nucleiChannel, nucleiBlurRadiusXY, nucleiBlurRadiusZ, nucleiMedian3DradiusXY, nucleiMedian3DradiusZ);
	else if (nucleiSegmentationChoice == "Cellpose cytoplasm segmentation") nuclei_info = segmentCytoplasmCellpose(process_image, nucleiChannel, probabilityThreshold, pixelWidth, unit, resultTable);
	labelmap_nuclei = nuclei_info[0];
	nuclei_outlines = nuclei_info[1];
	nrNuclei = nuclei_info[2];

	if(manualRemoveNuclei != "No thanks" && nrNuclei > 0) {
		nuclei_info = manually_remove_labels(labelmap_nuclei, nuclei_outlines, nrNuclei, process_image, imageName);
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
		nucleus_area_ = multiplyArraywithScalar(nucleus_area_, Math.sqr(pixelWidth));
		Table.setColumn("Nucleus ID", nucleus_id_, resultTable);
		Table.setColumn("Nucleus area 2D ("+unit+"^2)", nucleus_area_, resultTable);

		//Create a 3D version of the 2D nuclei labelmap
		if(gslices>1) Ext.CLIJ2_imageToStack(labelmap_nuclei, labelmap_nuclei_3D, gslices);
		else labelmap_nuclei_3D = labelmap_nuclei;

		//Foci filtering and detection - in a loop to enable parameter optimization 
		firstTimeProcessing = true;
	//	activeChannels = "111";
		zoom = 1;
		
		//Add selected ROI to the ROI Manager and add as overlay - must do this after StarDist
		if(useROI) {
			selectWindow(original);
			roiManager("reset");
			roiManager("add");
			roiManager("Select",0);
			roiManager("Set Color", "white");
			roiManager("Set Line Width", 1);
			roiManager("add");
		}
		
		do {
			//Check if fociChannelB actually exists and is different from Foci channel A
			if(detectFociChannelB == true && fociChannelB > gchannels) exit("Error: The selected Foci channel B ("+fociChannelB+") is higher than the number of channels of the image ("+gchannels+")."); 
			if(detectFociChannelB == true && fociChannelB == fociChannelA) {
				showMessageWithCancel("Selected Foci channels are the same", "Warning: The selected Foci channel B ("+fociChannelB+") is the same as Foci channel A!\nIf you continue only this channel will be analyzed.");
				detectFociChannelB = false;
			}

			//Foci channel A
			detections_fociA = detectFoci(process_image, fociChannelA, fociSizeA, anisotropyFactor, firstTimeProcessing, labelmap_nuclei, labelmap_nuclei_3D, thresholdFactorA);
			labelmap_fociA = detections_fociA[0];
			mask_fociA = detections_fociA[1];
			spots_fociA = detections_fociA[2];

			if(firstTimeProcessing == false) call("ij.gui.ImageWindow.setNextLocation", x_image, y_image);
			overlayA = mergeOriginalAndDetection(original, nrNuclei, nuclei_outlines, mask_fociA, spots_fociA, fociChannelA, useROI, x_ROI, y_ROI);
			overlay_image = overlayA;	//Will be overwritten if channel B is also used
			
			if(firstTimeProcessing == false && detectFociChannelB == false) close("Processing...");
			//Foci channel B
			if(detectFociChannelB == true) {
				detections_fociB = detectFoci(process_image, fociChannelB, fociSizeB, anisotropyFactor, firstTimeProcessing, labelmap_nuclei, labelmap_nuclei_3D, thresholdFactorB);
				labelmap_fociB = detections_fociB[0];
				mask_fociB = detections_fociB[1];
				spots_fociB = detections_fociB[2];
				overlayB = mergeOriginalAndDetection(original, nrNuclei, nuclei_outlines, mask_fociB, spots_fociB, fociChannelB, useROI, x_ROI, y_ROI);
				if(firstTimeProcessing == false) close("Processing...");
				Overlay.copy();	//Preserve nuclei outline overlays
				overlay_image = "Foci_overlay_ch"+fociChannelA+"_and_ch"+fociChannelB;
				run("Concatenate...", "  title="+overlay_image+" image1="+overlayA+" image2="+overlayB+" image3=[-- None --]");
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
//			Property.set("CompositeProjection", "Max");	//Use 'Composite Max' display setting
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
				Dialog.addChoice("Detection method", newArray("Marker-controlled watershed (recommended)","AreaMaxima local maximum detection"), foci_method);
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
				Dialog.addChoice("Nuclei outline visualization", newArray("bright","dim"), overlayBrightness);
				Dialog.addCheckbox("Done optimizing for this image only", false);
				Dialog.addCheckbox("Done optimizing; process all other images with these settings", false);

				//determine dialog location
				if(x_image + imageWidth + 500 < screenWidth) Dialog.setLocation(x_image+imageWidth-15, y_image);
				else if(x_image > 500) Dialog.setLocation(x_image-530, y_image);
//				else if(y_image + gheight + 300 < screenHeight) Dialog.setLocation(x_image, y_image+gheight+100);
				else Dialog.setLocation(x_image, y_image);
				Dialog.show();
//setBatchMode(true);
				//Get dialog entries
				fociSizeA = Dialog.getChoice();
				if(detectFociChannelB == true) fociSizeB = Dialog.getChoice();
				foci_method = Dialog.getChoice();
				thresholdFactorA = Dialog.getNumber();
				if(detectFociChannelB) thresholdFactorB = Dialog.getNumber();
				minFociSize = Dialog.getNumber();
				maxFociSize = Dialog.getNumber();
				if(detectFociChannelB) minOverlapSize = Dialog.getNumber();
				overlayBrightness = Dialog.getChoice();
				doneOptimizing = Dialog.getCheckbox();
				processAllOtherImages = Dialog.getCheckbox();
				if (processAllOtherImages == true) doneOptimizing = true;

				//get image display properties
				selectWindow(overlay_image);
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
					close(mask_fociA);
//					close("foci_RAW"+fociChannelA);
					if(detectFociChannelB == false) {
						selectWindow("Foci_overlay_ch"+fociChannelA);
						rename("Processing...");
					}
					else if(detectFociChannelB == true) {
						close("foci_ch"+fociChannelB);
						//close("MAX_Foci_ch"+channel+"_filtered_and_masked");
						close("Labelmap_detected_foci_filtered_ch"+fociChannelB);
						close(mask_fociB);
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

		if(useROI == false) {
			//Get nuclei positions
			selectWindow("Results");
			boundingBox_X = Table.getColumn("BOUNDING_BOX_X");
			boundingBox_Y = Table.getColumn("BOUNDING_BOX_Y");
			boundingBox_width = Table.getColumn("BOUNDING_BOX_WIDTH");
			boundingBox_height = Table.getColumn("BOUNDING_BOX_HEIGHT");
			Array.getStatistics(boundingBox_width, minWidth, maxWidth);
			Array.getStatistics(boundingBox_height, minHeight, maxHeight);
			if(debugMode) print("\nMaximum nucleus bounding box: "+maxWidth+", "+maxHeight);		
			
			//Draw nuclei numbers as overlay
			if(addNumbersOverlay) {
				setFont("SansSerif", labelFontSize, "bold antialiased");
				color = color_to_hex(fontColor);
				setColor(color);
				for (i = 0; i < nrNuclei; i++) {
					x = getResult("MASS_CENTER_X", i);
					y = getResult("MASS_CENTER_Y", i);
					Overlay.drawString(i+1, x - labelFontSize/2 + x_ROI, y + labelFontSize/2 + y_ROI);
				}
				updateDisplay();
			}
			//Measure the foci 
			measureFoci(original, fociChannelA, nrNuclei, labelmap_nuclei_3D, labelmap_fociA, boundingBox_X, boundingBox_Y, maxWidth, maxHeight);
			if(detectFociChannelB) measureFoci(original, fociChannelB, nrNuclei, labelmap_nuclei_3D, labelmap_fociB, boundingBox_X, boundingBox_Y, maxWidth, maxHeight);
	
	//		print("\n\nGPU Memory after channel "+fociChannelB);
	//		Ext.CLIJ2_reportMemory();
		
			//Compute A->B foci colocalization
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
			List.set("06. 3D handling?", ThreeDHandling);
			List.set("07. [single z-slice foci detection only] Slice nr", singleSlice);
			List.set("08. Remove image borders (pixels)", cropBorder);
			List.set("09. Image XY binning before analysis [0-n]", XYBinning);
			List.set("10. Nuclei segmentation method", nucleiSegmentationChoice);
			List.set("11. Nuclei binning factor [0-n] (StarDist)", downsampleFactorStarDist);
			List.set("12. Nuclei probablility threshold [0.0-1.0] (StarDist)", probabilityThreshold);
			List.set("13. Remove nulei with diameter smaller than (units)", minNucleusSize_setting);
			List.set("14. Remove nulei with diameter larger than (units)", maxNucleusSize_setting);
			List.set("15. Exclude nulei on image edges", excludeOnEdges);
			List.set("16. Enable foci parameters optimization mode?", optimizationMode);
			List.set("17. Foci size channel A", fociSizeA);
			List.set("18. Foci size channel B", fociSizeB);
			List.set("19. Foci detection method", foci_method);
			List.set("20. Foci intensity threshold bias channel A", thresholdFactorA);
			List.set("21. Foci intensity threshold bias channel B", thresholdFactorB);
			List.set("22. Minimum foci size", minFociSize);
			List.set("23. Maximum foci size", maxFociSize);
			List.set("25. Minimum overlap of foci to colocalize (pixels/voxels)", minOverlapSize);
			list = List.getList();
	
			Table.save(outputFolder + File.separator + imageName + "__Foci_results.tsv");
			if(saveOverlayImage == true) {
				if(detectFociChannelB == true) {
					selectWindow("Foci overlap map");
					run("RGB Color", "slices");
					selectWindow("Foci overlap map");
					Stack.setSlice(gslices/2);
					setBatchMode("show");
					setMetadata("info", list);
					saveAs("Gif", outputFolder + File.separator + imageName + "__foci_overlap_map");
				}
			}
			selectWindow(overlay_image);
			setMetadata("info", list);
			saveAs("tif", outputFolder + File.separator + imageName + "__" + overlay_image);
			run("Clear Results");

			Ext.CLIJ2_release(labelmap_nuclei);
			Ext.CLIJ2_release("foci_ch"+fociChannelA);
			if(detectFociChannelB == true) Ext.CLIJ2_release("foci_ch"+fociChannelB);
			Ext.CLIJ2_release("Labelmap_detected_foci_ch"+fociChannelA);
			if(detectFociChannelB == true) Ext.CLIJ2_release("Labelmap_detected_foci_ch"+fociChannelB);
	//		Ext.CLIJ2_reportMemory();
		}
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



function segmentNucleiClassic(image, channel, nucleiBlurRadiusXY, nucleiBlurRadiusZ, nucleiMedian3DradiusXY, nucleiMedian3DradiusZ) {
	selectWindow(image);
	getDimensions(width, height, channels, slices, frames);
	minNucleusSize = PI*Math.sqr((minNucleusSize_setting / pixelWidth / 2));	//Calculate the nucleus area as if it were a circle
	maxNucleusSize = PI*Math.sqr((maxNucleusSize_setting / pixelWidth / 2));
	
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

	Ext.CLIJ2_automaticThreshold(nuclei_filtered_MAX, thresholded, "Otsu");
	Ext.CLIJ2_pullBinary(thresholded);
	run("Watershed");	//Not on GPU becayse results are not so good.
	run("Properties...", "unit=&unit pixel_width=&pixelWidth pixel_height=&pixelHeight voxel_depth=1.0000");

	if(excludeOnEdges) run("Analyze Particles...", "size=" + minNucleusSize + "-" + maxNucleusSize + " circularity=0.25-1.00 show=[Count Masks] exclude include add");
	else run("Analyze Particles...", "size=" + minNucleusSize + "-" + maxNucleusSize + " pixel circularity=0.25-1.00 show=[Count Masks] include add");
	rename("Labelmap_nuclei");
	labelmap_nuclei = "Labelmap_nuclei";
	Ext.CLIJ2_push(labelmap_nuclei);
	Ext.CLIJ2_getMaximumOfAllPixels(labelmap_nuclei, nrNuclei);	//Count the number of nuclei

	//Create nuclei outlines from nuclei labelmap
	Ext.CLIJ2_detectLabelEdges(labelmap_nuclei, labelmap_edges);
	if(thickOutlines == false) {
		Ext.CLIJ2_mask(labelmap_nuclei, labelmap_edges, labelmap_outlines);
		Ext.CLIJ2_release(labelmap_edges);
	}
	else labelmap_outlines = labelmap_edges;
	Ext.CLIJ2_pullBinary(labelmap_outlines);
	Ext.CLIJ2_release(labelmap_outlines);
	rename("nuclei_outlines");
	run(outlineColor);
	
	nuclei_outlines = "nuclei_outlines";
	labelmap_nuclei = "Labelmap_nuclei";
	
	return newArray(labelmap_nuclei, nuclei_outlines, nrNuclei);
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
	getDimensions(dswidth, dsheight, dschannels, dsslices, dsframes);
	starDistTiles = pow(floor((maxOf(dswidth, dsheight)/maxTileSize)-1)+1,2);	//Determine the nr. of tiles
	run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':'"+starDist_input_image+"', 'modelChoice':'Versatile (fluorescent nuclei)', 'normalizeInput':'true', 'percentileBottom':'1.0', 'percentileTop':'99.60000000000001', 'probThresh':'"+probabilityThreshold+"', 'nmsThresh':'0.3', 'outputType':'ROI Manager', 'nTiles':'"+starDistTiles+"', 'excludeBoundary':'2', 'roiPosition':'Stack', 'verbose':'false', 'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");

	//Scale up ROIs
	if(downsampleFactorStarDist > 1) RoiManager.scale(downsampleFactorStarDist, downsampleFactorStarDist, false);

	//Spline fit ROIs
	selectWindow(image);	//Have to select an image with the correct size
	for (i = 0; i < roiManager("count"); i++) {
		roiManager("select", i);
		Roi.getSplineAnchors(x, y)
		Roi.setPolygonSplineAnchors(x, y);
		roiManager("Update");
	}
	
	//Convert ROIs to label map (using SCF function)
	run("Select None");
	run("ROI Manager to LabelMap(2D)");
	labelmap_nuclei = "Labelmap_nuclei_unfiltered";
	rename(labelmap_nuclei);
	run("Clear Results");

	selectWindow(labelmap_nuclei);
	getDimensions(uswidth, usheight, uschannels, usslices, usframes);
	if(uswidth != width || usheight != height) run("Canvas Size...", "width="+width+" height="+height+" position=Center zero");	// Make sure that the size of the upscaled labelmap is correct - upscaling can cause rounding errors

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
//	if(downsampleFactorStarDist > 1) {
//		labelmap_nuclei_final = "Labelmap_nuclei";
//		Ext.CLIJ2_greyscaleOpeningSphere(labelmap_nuclei_gapsClosed, labelmap_nuclei_final, downsampleFactorStarDist + 1, downsampleFactorStarDist + 1, 0);	//Smooth labels a bit
//		Ext.CLIJ2_release(labelmap_nuclei_gapsClosed);
//	}
//	else labelmap_nuclei_final = labelmap_nuclei_gapsClosed;
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
	run(outlineColor);
	
	nuclei_outlines = "nuclei_outlines";
	labelmap_nuclei_final = "Labelmap_nuclei";

	return newArray(labelmap_nuclei_final, nuclei_outlines, nrNuclei);
}



function segmentCytoplasmCellpose (image, channel, probabilityThreshold, pixelWidth, unit, resultTable) {
	selectWindow(image);
	
	getDimensions(width, height, channels, slices, frames);
	minNucleusSize = PI*Math.sqr((minNucleusSize_setting / pixelWidth / 2));	//Calculate the nucleus area as if it were a circle
	maxNucleusSize = PI*Math.sqr((maxNucleusSize_setting / pixelWidth / 2));

	run("Duplicate...", "title=nuclei duplicate channels=" +channel);	//Get nuclei (cells) for segmentation
	cellpose_input_image = getTitle();

	if(slices>1) {
		run("Z Project...", "projection=[Max Intensity]");
		rename("forCellpose_Zprojected");
		cellpose_input_image = getTitle();
	}

	setBatchMode("exit and display"); //Cellpose doesn't work in batch mode
	selectWindow(cellpose_input_image);
	run("Cellpose Advanced", "diameter=0 cellproba_threshold=0 flow_threshold="+probabilityThreshold+" anisotropy=1.0 diam_threshold=12.0 model=cyto nuclei_channel=0 cyto_channel=1 dimensionmode=2D stitch_threshold=-1.0 omni=false cluster=false additionnal_flags=-");
	labelmap_nuclei = getTitle();
//	selectWindow(image);
	selectWindow(labelmap_nuclei);
	setBatchMode("hide");
	selectWindow("nuclei");
	setBatchMode("hide");
	
	//Close unused images
	if(slices>1) close("forCellpose_Zprojected");

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
	run(outlineColor);
	
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
	setBatchMode("show");

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
	else if(debugMode) showImagefromGPU(foci);

	//Determine filter sizes
	if(fociSize == "tiny")	{ fociFilterSizeXY = 0;	fociSizeXY = 0.5;		}
	if(fociSize == "small")	{ fociFilterSizeXY = 0.5;	fociSizeXY = 1;		}
	if(fociSize == "average"){ fociFilterSizeXY = 1.0;	fociSizeXY = 1.5;	}
	if(fociSize == "large")	{ fociFilterSizeXY = 2;	fociSizeXY = 2;			}
	if(fociSize == "huge")	{ fociFilterSizeXY = 4;	fociSizeXY = 4;			}
	if(fociSize == "other") {
		fociSize = getNumber("Enter foci size in pixels", 2);
		fociFilterSizeXY = fociSize;
		fociSizeXY = fociSize;
	}

	//pull from GPU, remove outliers, push to GPU
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
	print("Detecting foci in channel "+channel); 
	print("Median of the Standard Deviations of all nuclei: "+medianNucleiStddev);

	fociFilterSizeZ = 2*anisotropyFactor * fociFilterSizeXY;	//semi-arbitrary setting
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
		//if(bits == 8) run("16-bit");
	}
//TO DO: Check Labelmap_detected_foci for 8-bit images? (no values >255)
	threshold = medianNucleiStddev * thresholdMultiplier;	//Set default threshold at n times the stddev (of the outlier-removed foci, which is a bit lower than the actual stddev, but less dependent on many foci being present)
	threshold = threshold*exp(thresholdFactor);
	print("Threshold set at "+d2s(thresholdMultiplier*exp(thresholdFactor),1)+" times above background stddev (bias set at "+thresholdFactor+"): Threshold used: "+threshold);
	maxDisplaySetting = minOf(pow(2,bits), threshold * 5);
//logWindow = split(getInfo("log"),"\n");
//waitForUser(logWindow.length);
	if(foci_method == "AreaMaxima local maximum detection") {
		Ext.CLIJ2_pull(foci_filtered);
		selectWindow(foci_filtered);
//		setBatchMode("show");
//		hMin = 100;
		//setBatchMode(false);	//H-watershed does not work in batch mode if this image is not shown.
//		run("H_Watershed", "impin=["+foci+"] hmin="+hMin+" thresh="+threshold+" peakflooding="+H_max_peak_flooding + " outputmask=false allowsplitting=false");
//		print(getTitle);
		run("AreaMaxima local maximum detection (2D, 3D)", "minimum="+minFociSize+" threshold="+threshold/2);
		rename("Labelmap_detected_foci_ch"+channel);
		Ext.CLIJ2_push("Labelmap_detected_foci_ch"+channel);

		labeledSpots = "labeledSpots_ch"+channel;
		Ext.CLIJ2_reduceLabelsToCentroids("Labelmap_detected_foci_ch"+channel, labeledSpots);
		Ext.CLIJ2_mask(labeledSpots, labelmap_nuclei_3D, labeledSpots_masked);
		if(isOpen("foci_spots_ch"+channel)) close("foci_spots_ch"+channel);
		Ext.CLIJ2_pullBinary(labeledSpots_masked);
		rename("foci_spots_ch"+channel);
		Ext.CLIJ2_release(labeledSpots);
		Ext.CLIJ2_release(labeledSpots_masked);
		close(foci_filtered);
		close("Labelmap_detected_foci_ch"+channel);
	}

	else if(foci_method == "Marker-controlled watershed (recommended)") {
		//	Alternative function - could be faster: 
		//	Ext.CLIJx_detectAndLabelMaximaAboveThreshold(foci, detectedMaxima, 0, 0, 0, threshold, false)
		//	Ext.CLIJ2_pull(detectedMaxima);
		//	rename("detectedMaxima");

		//Detect all maxima
		allMaxima = "allMaxima";
		if(ThreeDHandling == "Use quasi-2D foci detection (detect foci separately in every Z-slice)" && gslices>1) Ext.CLIJ2_detectMaximaSliceBySliceBox(foci_filtered, allMaxima, fociSizeXY, fociSizeXY);
		else Ext.CLIJ2_detectMaxima3DBox(foci_filtered, allMaxima, fociSizeXY, fociSizeXY, fociSizeZ);
		//remove spots lower than threshold
		MaskFociAboveThreshold = "MaskFociAboveThreshold_ch"+channel;
		Ext.CLIJ2_threshold(foci_filtered, MaskFociAboveThreshold, threshold);
		Ext.CLIJ2_mask(allMaxima, MaskFociAboveThreshold, maskedSpots);
		Ext.CLIJ2_release(allMaxima);
		if(debugMode) showImagefromGPU(MaskFociAboveThreshold);
		labeledSpots = "labeledSpots_ch"+channel;
		Ext.CLIJ2_connectedComponentsLabelingBox(maskedSpots, labeledSpots);	//Create label image with single pixels
		Ext.CLIJ2_mask(labeledSpots, labelmap_nuclei_3D, labeledSpots_masked);
		Ext.CLIJ2_release(maskedSpots);
		Ext.CLIJ2_release(labeledSpots);
//		if(bits == 8) Ext.CLIJ2_release(maskedSpots_16bit);
		if(debugMode) {
			showImagefromGPU(labeledSpots_masked);
			//run("16-bit");
			run("glasbey on dark");
			run("Merge Channels...", "c1="+foci_filtered+" c2="+labeledSpots_masked+" create keep");
			rename("Overlay_spots_ch"+channel);
			Stack.setChannel(2);
			run("Magenta");
			Stack.setChannel(1);
			run("Green");
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
		Ext.CLIJx_morphoLibJMarkerControlledWatershed(foci_filtered_inverted, labeledSpots_masked, MaskFociAboveThreshold, labelmap_foci);
		Ext.CLIJx_morphoLibJFillHoles(labelmap_foci, labelmap_foci_filled);
		Ext.CLIJ2_release(labelmap_foci);

		if(isOpen("foci_spots_ch"+channel)) close("foci_spots_ch"+channel);
		Ext.CLIJ2_pullBinary(labeledSpots_masked);
		rename("foci_spots_ch"+channel);
		Ext.CLIJ2_release(foci_filtered);
		Ext.CLIJ2_release(foci_filtered_inverted);
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

	//Remove some lines in the log window (unnecessary info by MorphoLibJ)
	logWindowContents = getInfo("log");
	logWindowLines = split(logWindowContents, "\n");
	logWindowLines = Array.trim(logWindowLines, logWindowLines.length - 7);
	logWindowLines = arrayToString(logWindowLines);
	print("\\Clear");
	print(logWindowLines);

	print("\\Update:Channel "+channel+": "+max+" foci detected in "+nrNuclei+" nuclei ("+d2s(max/nrNuclei,1)+" foci per nucleus)" + "\n");

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
	if(isOpen("foci_spots_ch"+channel)) close("Mask_foci_ch"+channel);
	Ext.CLIJ2_pullBinary(labelmap_foci_final);
	rename("Mask_foci_ch"+channel);

	Ext.CLIJ2_release(labelmap_foci_final);
	if(!debugMode) close(foci+"_outliers_removed");

	return newArray("Labelmap_detected_foci_filtered_ch"+channel, "Mask_foci_ch"+channel, "foci_spots_ch"+channel);
}


function mergeOriginalAndDetection(original, nrNuclei, nuclei_outlines, foci_mask, foci_spots, channel, useROI, x_ROI, y_ROI) {
	// Merge the original foci channel with the detected foci

	if(useROI == false) {
		//prepare nuclei
		selectWindow("nuclei");
		run(bits+"-bit");	//convert to allow merging
	
		//prepare foci raw
		foci_RAW = "foci_ch"+channel;
		if(bits == 16) Ext.CLIJ2_convertUInt16(foci_RAW, foci_RAW_8or16bit);
		else if(bits == 8) Ext.CLIJ2_convertUInt8(foci_RAW, foci_RAW_8or16bit);
		Ext.CLIJ2_pull(foci_RAW_8or16bit);
		Ext.CLIJ2_release(foci_RAW_8or16bit);
	}
	if(useROI == true) {
		//prepare nuclei - channel 1
		close("nuclei");
		selectWindow(original);
		run("Select None");
		run("Duplicate...", "title=nuclei duplicate channels="+nucleiChannel);
		run(bits+"-bit");	//convert to allow merging
		
		//prepare foci raw - channel 2
		selectWindow(original);
		foci_RAW_8or16bit = "foci_ch"+channel;
		run("Duplicate...", "title="+foci_RAW_8or16bit+ " duplicate channels="+channel);
		run(bits+"-bit");	//convert to allow merging
	}
	
	//prepare foci_mask - channel 3
	selectWindow(foci_mask);
	run(bits+"-bit");	//convert to allow merging
	if(useROI == true) {
		Ext.CLIJ2_push(foci_mask);
		Ext.CLIJ2_create3D(foci_mask_large, gwidth, gheight, gslices, bits);
		Ext.CLIJ2_set(foci_mask_large, 0);
		Ext.CLIJ2_paste3D(foci_mask, foci_mask_large, x_ROI, y_ROI, 0);
		Ext.CLIJ2_pull(foci_mask_large);
		setMinAndMax(0, 255);
	}
	else foci_mask_large = foci_mask;
	
	//prepare foci_spots - channel 4
	selectWindow(foci_spots);
	run(bits+"-bit");	//convert to allow merging
	if(useROI == true) {
		Ext.CLIJ2_push(foci_spots);
		Ext.CLIJ2_create3D(foci_spots_large, gwidth, gheight, gslices, bits);
		Ext.CLIJ2_set(foci_spots_large, 0);
		Ext.CLIJ2_paste3D(foci_spots, foci_spots_large, x_ROI, y_ROI, 0);
		Ext.CLIJ2_pull(foci_spots_large);
		setMinAndMax(0, 255);
	}
	else foci_spots_large = foci_spots;
	
if(debugMode) {
	selectWindow("nuclei");
	setBatchMode("show");
//	selectWindow(foci_RAW_8or16bit);
//	setBatchMode("show");
	selectWindow(foci_mask_large);
	setBatchMode("show");
	selectWindow(foci_spots_large);
	setBatchMode("show");
}
//setBatchMode("exit and display");
//exit;
	run("Merge Channels...", "c1=nuclei c2=" + foci_RAW_8or16bit + " c3=" + foci_mask_large + " c4=" + foci_spots_large + " create keep");
	rename("Foci_overlay_ch"+channel);
	getDimensions(width, height, channels, slices, frames);
	Stack.setSlice(slices/2);
	Stack.setChannel(1);
	setLabel("Foci_overlay_ch"+channel, "NUCLEI");
    setLut(reds, greens, blues);
 	run("Enhance Contrast...", "saturated=0.02");  
	Stack.setChannel(2);
	setLabel("Foci_overlay_ch"+channel, "FOCI");
	run("Green");
	run("Enhance Contrast...", "saturated=0.02");
	getMinAndMax(min, max);
	setMinAndMax(min, maxOf(max, min + maxDisplaySetting));
	Stack.setChannel(3);
	setLabel("Foci_overlay_ch"+channel, "DETECTED FOCI");
	resetMinAndMax;
	run("Magenta");
	Stack.setChannel(4);
	setLabel("Foci_overlay_ch"+channel, "CENTROIDS OF DETECTED FOCI");
	resetMinAndMax;
	run("Grays");	

	//Add nuclei outlines as overlay to the original image
	selectWindow(nuclei_outlines);
	if(overlayBrightness == "dim") setMinAndMax(0, 512);	//Silly way to set nuclei outline overlay to dim, but it works
	if(overlayBrightness == "bright") resetMinAndMax;	//Silly way to set nuclei outline overlay to dim, but it works
	if(debugMode) setBatchMode("show");
	selectWindow("Foci_overlay_ch"+channel);
	for (i = 1; i <= gslices; i++) {
		if(gslices>1) Stack.setSlice(i);
		run("Add Image...", "image="+nuclei_outlines+" x=0 y=0 opacity="+labelOpacity+" zero");	//Add labelmap to image as overlay
		Overlay.setPosition(0, 0, 0);
	}
/*
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
*/
	if(useROI) {
		run("Translate...", "x="+x_ROI+" y="+y_ROI+" interpolation=None overlay");
		run("From ROI Manager");	//Add the selected ROI as overlay
	}

	return "Foci_overlay_ch"+channel;
}


function setLabel(image, label) {
	selectWindow(image);
	getDimensions(width, height, channels, slices, frames);
	Stack.getPosition(ch, sl, fr);
	if(slices>1) {
		for (i=1; i<=slices; i++) {
			Stack.setSlice(i);
			//run("Set Label...", "label="+label.replace(" ","_"));
			setMetadata("Label", label.replace(" ","_"));
		}
	}
	else setMetadata("Label", label.replace(" ","_"));
	Stack.setPosition(ch, sl, fr);
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
		
		if(i%10 == 0) run("Clear Results");	//A large Results table is slowing down a lot! Intermediate clearing helps.
	}
	Ext.CLIJ2_release(foci_raw_cropped);
	Ext.CLIJ2_release(labelmap_nuclei_cropped);
	Ext.CLIJ2_release(labelmap_foci_cropped);
	Ext.CLIJ2_release(labelmap_foci_nucleus_relabeled);
	Ext.CLIJ2_release(labelmap_foci_nucleus);
	if(debugMode) print("Measuring all foci in channel "+channel+" took "+ getTime() - startTime +" ms.");

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
		if(c == fociChannelA || detectFociChannelB == true && c == fociChannelB) continue;
		else {
			selectWindow(original);
			run("Clear Results");
			Stack.setChannel(c);
			channel_to_measure = original+"_ch"+c;
			run("Duplicate...", "title=[" + channel_to_measure + "] duplicate channels=" + c);
			Ext.CLIJ2_push(channel_to_measure);
			Ext.CLIJ2_statisticsOfBackgroundAndLabelledPixels(channel_to_measure, labelmap_nuclei_3D);
			mean_intensity_ = Table.getColumn("MEAN_INTENSITY", "Results");
			sum_intensity_ = Table.getColumn("SUM_INTENSITY", "Results");
			background = mean_intensity_[0];
			Array.fill(background_, background);
			mean_intensity_ = Array.slice(mean_intensity_, 1, mean_intensity_.length);
			sum_intensity_ = Array.slice(sum_intensity_, 1, sum_intensity_.length);

			Table.setColumn("Background intensity ch"+c, background_, resultTable);
			Table.setColumn("Mean nucleus intensity ch"+c, mean_intensity_, resultTable);
			Table.setColumn("Sum nucleus intensity ch"+c, sum_intensity_, resultTable);
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

	if(debugMode) print("Measuring colocalization took "+ getTime() - startTime +" ms.");

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

//Create a single color LUT for the nuclei overlay
function create_lut(reds, greens, blues) {
	for (i=0; i<256; i++) {
		reds[i] = i;
		greens[i] = i;
		blues[i] = i;
	}
	return newArray(reds, greens, blues);
}

//Convert a color into a hexadecimal code
function color_to_hex(color) {
	colorArray = split(color,",,");
	hexcolor = "#" + IJ.pad(toHex(colorArray[0]),2) + IJ.pad(toHex(colorArray[1]),2) + IJ.pad(toHex(colorArray[2]),2);
	return hexcolor;
}


//Converts an array into a string, elements separated by "\n"
function arrayToString(array) {
	outputString = "";
	for (i = 0; i < array.length; i++) {
		outputString += array[i] + "\n";
	}
	return outputString;
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

