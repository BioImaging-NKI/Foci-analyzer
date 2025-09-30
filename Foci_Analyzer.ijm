/* Macro to quantify foci in nuclei/cells. Works on 2D/3D images, including multiseries files, as long as all series have the same specs
 * More info on https://imagej.net/plugins/foci-analyzer
 * 
 * ► Requires the following Fiji update sites:
 * - 3D ImageJ Suite
 * - CLIJ
 * - CLIJ2
 * - CLIJx-assistent
 * - CLIJx-assistent-extensions
 * - CSBDeep
 * - IJPB-plugins
 * - StarDist
 *(- TensorFlow) In case StarDist gives an error. See https://forum.image.sc/t/stardist-error-since-update/107729
 *
 * ► In order to run Cellpose segmentation you also need:
 * - A working Cellpose Python environment
 * - PTBIOP update site, with proper settings. See https://github.com/BIOP/ijl-utilities-wrappers/blob/master/README.md#cellpose
 * 
 * 
 * Author: Bram van den Broek, The Netherlands Cancer Institute
 * For questions please use the Image.sc forum (https://forum.image.sc/) with tag @bramvdbroek.
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
 * Version 1.00, July 2022:
 * - Added possibility to manually remove segmented nuclei by clicking. Masks are saved and can be loaded for re-analysis.
 * 
 * Version 1.10, July 2022:
 * - Fixed a small bug: nuclear intensity in foci channel B is now measured also when detectFociChannelB is false.
 * - Changes in AreaMaximum detection (minimum foci size) 
 * - Now asks to create output folder if it doesn't exist.
 * - Classic segmentation works again
 * - Fixed a bug where the nucleus area was not correctly reported when Stardist nuclei binning factor > 1
 * 
 * Version 1.20, October 2022:
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
 * Version 1.30, December 2022
 * - If 'Alt' is pressed before a file is loaded, the user can select a ROI to do a preview analysis on.
 * - downscaled ROIs are now smoothed using spline fitting
 * - Improved GUI: 3D image handling radiobuttons instead of separate checkboxes
 * - Improved logging: Threshold factor * bias is now mentioned in the log window
 * - Colocalization overlap image is now also 3D and saved as .Gif file 
 * 
 * Version 1.32, April 2023
 * - Added the actually used threshold value to the output image metadata
 * - Included option to detect foci on extended depth of focus projection (Sobel, GPU)
 * 
 * Version 1.35, June 2023
 * - Added option to include foci outside nuclei. Distance range can be specified.
 * 
 * Version 1.36, February 2024
 * - Fixed a bug causing a crash with 32-bit images when producing the overlay.
 * - Added a cellpose diameter parameter
 * 
 * Renamed to version 1.4 for release on GitHub: https://github.com/BioImaging-NKI/Foci-analyzer/
 *
 * Version 1.41, May 2024
 * - Added possibility to load 2D nuclei segmentations from a .zip file containing ImageJ ROIs (e.g. exported from QuPath)
 * 
 * Version 1.5, September 2024
 * - Added foci overlap count and area/volume for both channels to the result, as well as in the overlap map image (thanks to Mabel Baxter Dalrymple) 
 * - Overlay and overlap images are now saves as .zip files to save storage space (typically 10 times), at the expense of slower processing.
 * - Possibility for automatic rescaling nuclei for StarDist (to a pixel size 0.4 µm), with non-integer rescaling factor.
 *   Upscaling (downscale factor < 1) is now also allowed. Interpolation is set to bilinear.
 * - Settings optimization: Channels Tool and B&C window are displayed in upper left corner of the screen. Added some tips in the optimization dialog.
 * 
 * Version 1.52, September 2024
 * - Added a 'Help' button in the optimization dialog, linking to the GitHub site.
 * - Added the possibility (in the optimization dialog) to use the same fixed absolute threshold for all images.
 * - Improved starting dialog, with headings
 * 
 * Version 1.53, November 2024
 * - Added possibility to fix the absolute threshold value for all images. Changed the optimization dialog for this.
 * - Channels and slices can now be flipped in images with multiple slices but only a single channel.
 * 
 * Version 1.54, December 2024
 * - Expanded Cellpose detection options.
 * 
 * Version 1.59, December 2024
 * - Measure foci using MorphoLibJ instead of CLIJ2. Faster and a bit easier to handle.
 * - A table with all foci statistics is now also saved (thanks to Harry Osborne)
 * - Added possibility to load settings from a previously analyzed file
 * 
 * Version 1.60-1.70, January-March 2025
 * - Cellpose envPath and envType parameters are automatically retrieved (new Cellpose wrapper)
 * - Added 3D segmentation with Cellpose 2.5D/3D, and 3D label visualization. Cellpose 3.1.0 is required (for --dP_smooth/flow3D_smooth parameter)
 * - Changed automatic StarDist downscaling to 0.25 um/pixel (was 0.4)
 * - Persistent foci detection script parameters are now saved after running the optimization dialog
 * - Fixed a bug where optimization dialog would be skipped (thanks to Harry Osborne)
 * - Foci centroids of removed foci due to size restrictions are now removed as well
 * - Improved filtering in Z -> better foci detection in 3D
 * - Fixed bug where the median foci volume was given in pixels, thanks to Harry Osborne   
 * - Fixed bug where the foci area was incorrect for 3D images when analyzed in 2D (foci area was multiplied by voxelDepth)
 * - Input images can now be time-lapse series. Frames are split and processed sequentially.
 * - Option to load label images instead of performing segmentation. (Label images need to have the same name as the original images, in a different folder.)
 * - If ROI/label images is not found it skips the image.
 * - Foci count can be overlayed on the image, as well as the nucleus number
 * - Added StdDev projection method for handling 3D images
 * - Use own code for ROI manager to labelmap (instead of SCF MPI CBG plugin)
 * 
 * Version 1.74:
 * - Possibility to choose a color for foci count overlay
 * - Possibility to choose 'none' for nucleus ID / foci count overlay
 * 
 * Version 1.80:
 * - Added checkbox to hide the 3D segmentation dialog window for subsequent images.
 * - Fixed bug where the macro would crash when 3D Cellpose segmentation was chosen in combination with *not* excluding cells on edges.
 * - Release version for Elmi 2025
 * 
 * Version 1.83:
 * - Added foci centroid coordinates to the 'All foci statistics' table (in 1.82 - bugfix in 1.83)
 * - Added descriptions to the Scijava script parameters
 * 
 * Version 1.84:
 * - Changed comments
 * 
 * Version 1.85:
 * - Fixed a bug that caused foci detection outside nuclei to be inaccurate (dilated labels were incorrectly cropped)
 * 
 * Version 1.86-1.87:
 * - New option to set file extension of images to be processed
 * - New option to save the output files in the same (sub)folder as the input files.
 * - Only perform StarDist and CSBDeep plugin check when using StarDist is actually used
 * - Small bugfixes regarding manual removal of nuclei, circumventing the issue that sometimes old multipoint selections stay in memory (v1.87)
 * 
 * Version 1.88:
 * - Implemented a 'Manual threshold' method for foci detection
 * - Gaussian filter in classic segmentation is now pixelsize-dependent (0.5 µm, or 2 pixels if unit is not microns)
 * - Fixed bug in classic segmentation where the nucleus min and max sizes were incorrect for calibrated images
 * - Updated Cellpose --dPSmooth parameter to --flow3D_smooth
 * 
 */

#@ String	Foci_Analyzer_message 	(value="<html><p style='font-size:18px; color:#000000; font-weight:bold'><img width=96 height=96 src='https://imagej.net/media/icons/Foci-Analyzer-icon.png'</img><a style='color:#000000' href=https://imagej.net/plugins/foci-analyzer>Foci Analyzer</a> (v1.88)</p></html>", visibility="MESSAGE")
#@ String	file_message 			(value="<html><p style='font-size:14px; color:#9933cc; font-weight:bold'>File settings</p></html>", visibility="MESSAGE")
#@ File[]	files 					(label = "Input files", style="File", description="Here you can specify which files to analyze, by adding them to the list, or drag&drop from a file explorer window.")
#@ String	processOnlyExtension	(label = "Only process files with extension (leave empty for all files)", value="", description="If files with multiple formats are in the list, only files ending with this extension (e.g. tif, czi) will be processed.")
#@ File		outputFolder			(label = "Output folder", style = "directory", required=false, description="The folder where all the analyzed images and results will be saved.")
#@ Boolean	useInputAsOutput		(label = "Save the output files in the same folder as the input file(s)?", description="When checked, output images and result files will be saved in the same folder as the input file(s).")

#@ String	image_message 			(value="<html><p style='font-size:14px; color:#3366cc; font-weight:bold'>Image settings</p></html>", visibility="MESSAGE")
#@ Boolean	loadSettingsFromFile	(label = "Load settings from previously analyzed image?", description="When checked, a separate dialog will popup where an output .zip file can be selected (overlay or colocalization map).\nAll settings are loaded from the metadata in the saved image. The script parameter entries in this large dialog are ignored.")
#@ Integer 	nucleiChannel			(label = "Nuclei channel (-1 if not used)", value = 1, description="The image channel that contains the nuclei. For StarDist nuclei segmentation is performed using this channel, which always happens in 2D (in the case of 3D images on a maximum intensity projection).\nFor Cellpose, there are two possibilities, depending on the value of Cytoplasm/membrane channel below.")
#@ Integer	cytoChannel 			(label = "Cytoplasm/membrane channel (-1 if not used)", value = -1, description="The image channel that contains cytoplasm or membrane. If set to -1, segmentation is performed on the nucleus channel alone.\nIf not -1 and Cellpose is chosen as segmentation method, segmentation is performed on this channel. In this case the nuclei channel is the 'additional helper channel'.\nIf not set to -1 and StarDist is chosen, this channel is not used in the segmentation, but instead just displayed in the overlay image. (default: -1)")
#@ Integer	fociChannelA 			(label = "Foci channel A", value = 2, description="The first foci channel (default: 2)")
#@ Integer	fociChannelB 			(label = "Foci channel B", value = 3, description="The second foci channel (default: 3)")
#@ Boolean	detectFociChannelB		(label = "Also detect foci in channel B and perform colocalization?", persistence=true, value=true, description="If checked, foci in both channels A and B will be analyzed, followed by a simple colocalization analysis. (default: checked)")

#@ String	ThreeDHandling 			(label = "2D/3D foci detection", choices={"Detect foci in 3D (or 2D when N/A)", "Detect foci on the Maximum Intensity Projection", "Detect foci on a Extended Depth of Focus Projection", "Detect foci on the Standard Deviation Projection", "Detect foci on the Summed Intensity Projection", "Use quasi-2D foci detection (detect foci separately in every Z-slice)", "Process a single z-slice only (specify below which slice)"}, style="listbox", description="2D/3D foci detection:\nThis parameter determines how foci in 3D images should be analyzed. For 2D input images this setting is ignored. The options are:\n\n- Detect foci in 3D (or 2D when N/A) (default) performs foci analysis using 3D filters and 3D marker-controlled watershed functions. Connected foci in consecutive slices are counted once.\n\n- Detect foci on the Maximum Intensity Projection performs 2D foci analysis on the MIP projection.\n\n- Detect foci on a Extended Depth of Focus Projection performs 2D foci analysis on an EDF projection.\n\n- Detect foci on the Standard Deviation Projection performs 2D foci analysis on the STDEV projection.\n\n- Detect foci on the Summed Intensity Projection performs 2D foci analysis on the SUM projection.\n\n- Use quasi-2D foci detection (detect foci separately in every Z-slice) analyzes every z-slice in a 3D image as a separate 2D image.\n  This setting is useful in cases where the z-spacing is very large and each focus is visible in only one z-slice.\n  Hence, connected foci in consecutive slices will be counted multiple times.\n\n- Process a single z-slice only (specify below which slice) allows the user to analyze foci only in a particular z-slice.")
#@ Integer	singleSlice				(label = "[single z-slice foci detection only] Slice nr", value=1, min=1, description="The single z-slice used for the option above. For any other choice this parameter is ignored.")
#@ Integer	cropBorder				(label = "Remove image borders (pixels)", value = 0, min=0, description="Possibility to remove edges from the image. This can in particular be useful when the image edges have very different intensities, causing incorrect automatic nuclei segmentation. (default: 0)")
#@ Integer	XYBinning				(label = "Image XY binning before analysis [1-n]", value = 1, min=1, description="Optional pixel binning in case the resolution is very high and the foci consist of many pixels.\nA value of 2 means: 2x2 pixels will be binned into 1 pixel. This reduces noise in the image and speeds up analysis. (default: 1)")

#@ String	nuclei_message 			(value = "<html><p style='font-size:14px; color:#33aa00; font-weight:bold'>Nuclei/cell detection settings</p></html>", visibility="MESSAGE")
#@ String	nucleiSegmentationChoice	(label = "Nuclei/cell segmentation method", choices={"StarDist nuclei segmentation 2D (or on 2D projection)", "Cellpose segmentation 2D (or on 2D projection)", "Cellpose segmentation 3D", "Classic segmentation", "Load ROIs from file", "Load label images"}, style="listBox", description="Nuclei/cell segmentation method:\nStardist nuclei segmentation 2D (or on 2D projection) (default) uses the pretrained convolutional neural network StarDist to recognize cell nuclei in fluorescence microscopy images.\nIn general this method works very well on a large variety of samples.\n\n- Cellpose segmentation 2D (or on 2D projection) uses the deep learning network Cellpose to recognize whole cells or nuclei.\nUse this option if you want to measure foci in entire cells, or if you prefer Cellpose nuclei segmentation over StarDist.\nN.B. Cellpose requires additional installations (see Installation / Requirements).\n\n- Cellpose segmentation 3D: If this option is chosen a new dialog pops up with extra settings. These are the most important parameters for 3D segmentation.\nMore parameters can be added in the 'Additional Cellpose parameters' field. The Help button takes you to the Cellpose CLI with explanations of all parameters.\n\n- Classic nuclei segmentation allows the user to segment nuclei using manual/automatic thresholding is provided for legacy reasons.\nThe method is almost always outperformed by the two other methods.\n\n- Load ROIs from file: ImageJ ROI .zip files can be loaded instead of performing segmentation. This option is used in the (near future) QuPath-Fiji workflow.\nROI files should have the same name as the input images without extensions, followed by '_ROIs.zip'.\n\n- Load label images allows loading a labelmap, if the segmentation has been done by external programs, or to quickly re-run files with the same segmentation.\nLabel image files should be present in another folder and have the exact same name as the input images.")
#@ Double	downsampleFactorStarDist	(label = "Stardist nuclei rescaling factor [1-n], 0 for automatic, 1 for no rescaling", value = 0, min=0, description="Stardist is trained on medium resolution images, and generally performs well on images with pixel sizes around 0.2-0.5 µm.\nSet to 0 for automatic rescaling the nuclei to an optimal pixel size of 0.25 µm, or put any other number for manual control of the rescaling.")
#@ Double 	probabilityThreshold		(label = "Probability/flow threshold [0.0-1.0] (StarDist/Cellpose)", value = 0.5, min=0, max=1, style="format:0.0", description="Lower values will accept more nuclei/cells; higher values will be more stringent. For Cellpose this is actually the flow_threshold parameter.")
#@ String 	CellposeModel			(label = "Cellpose model", choices={"cyto3","nuclei","tissuenet_cp3","cpsam","custom"}, style="listBox", value="cyto3", description="The model (built-in or custom) used for segmentation.")
#@ Integer	CellposeDiameter		(label = "Cellpose cell diameter (pixels), 0 for automatic", value = 0, min=0, description="Estimated diameter of the cells, in pixels. Setting this parameter to 0 will trigger Cellpose to estimate it.")

#@ Integer 	minNucleusSize_setting	(label = "Remove nuclei/cells with diameter smaller than (µm)", value = 4, description="Objects smaller than circles having an area corresponding to this diameter will be removed.")
#@ Integer	maxNucleusSize_setting	(label = "Remove nuclei/cells with diameter larger than (µm)", value = 50, description="Likewise, but for large objects.")
#@ Boolean	excludeOnEdges			(label = "Exclude nuclei/cells on image edges", value = true, description="When checked, nuclei that touch the image edge will not be analyzed. (default: checked).")
#@ String	manualRemoveNuclei		(label = "Manually remove segmented nuclei/cells", choices={"No thanks","Manually remove nuclei", "Load previously saved removals (from output folder)"}, value = "No thanks", description="Manual nuclei removal: allows the user to erase ill-segmented nuclei before analysis. (default: No thanks)\nOptions:\n\n- No thanks means no manual nuclei editing\n\n- Manually remove nuclei : Remove nuclei by leftclicking them in the image with the mouse. Editings will be saved to a small text file in the output folder.\n\n- Load previously saved removals (from output folder) : If you have edited the segmented nuclei on this image before, it will load the previous nuclei removals\n  from the file in the specified output folder. (Hence, if you change the output folder parameter this option will not work.)")

#@ String	foci_message1			(value="<html><p style='font-size:14px; color:#cc9933; font-weight:bold'>Foci detection settings</p></html>", visibility="MESSAGE")
#@ Boolean	optimizationModeSetting	(label = "Preview foci detection (for parameter fine-tuning)?", value=true, description="Checking this will allow the user to adapt the foci detection settings on a preview analysis before quantifying.")
#@ String	fociSizeA				(label = "Foci size channel A (after XY binning)", choices={"tiny (0.5 pixels)","small (1 pixel)","average (2 pixels)","large (4 pixels)","huge (8 pixels)","other (define later)"}, style="listBox", value="average", description="This parameter controls several foci image filtering steps and steers the macro towards detecting smaller or larger foci.")
#@ String	fociSizeB				(label = "Foci size channel B (after XY binning)", choices={"tiny (0.5 pixels)","small (1 pixel)","average (2 pixels)","large (4 pixels)","huge (8 pixels)","other (define later)"}, style="listBox", value="average", description="This parameter controls several foci image filtering steps and steers the macro towards detecting smaller or larger foci.")
#@ String	fociDetectionMethod		(label = "Foci detection method", choices={"Marker-controlled watershed (recommended)","AreaMaxima local maximum detection","Manual thresholding"}, style="listBox", description="Select the method for foci detection.")

#@ Double	thresholdFactorA		(label = "Foci intensity threshold bias channel A", value = 0.0, min=-2.5, max=2.5, style="scroll bar", stepSize=0.01, description="The macro will automatically estimate the intensity threshold for foci detection. This default threshold can be biased with the slider.")
#@ Double	thresholdFactorB		(label = "Foci intensity threshold bias channel B", value = 0.0, min=-2.5, max=2.5, style="scroll bar", stepSize=0.01, description="The macro will automatically estimate the intensity threshold for foci detection. This default threshold can be biased with the slider.")
#@ Integer	minFociSize				(label = "Minimum foci size (area) (pixels/voxels)", value = 3, description="Foci occupying an area/volume smaller than this value will be deleted.")
#@ Integer	maxFociSize				(label = "Maximum foci size (area) (pixels/voxels)", value = 9999, description="The upper limit for the foci size, in pixels/voxels.")

#@ Double	maxFociDistanceOutsideNuclei_setting	(label = "Max distance of foci outside nuclei/cells (µm); -1 for full image", value = 0, min=-1, style="format:0.0", description="This controls how far outside the cell/nucleus segmentation foci should still be counted. (default: 0)")

#@ Integer	minOverlapSize			(label = "Minimum overlap of foci to colocalize (pixels/voxels)", min = 1, value = 1, description="Foci in channel A and B will be counted as colocalizing only if they overlap with at least this area/volume.")

#@ String	visualization_message	(value = "<html><p style='font-size:14px; color:#cc3333; font-weight:bold'>Visualization options</p></html>", visibility="MESSAGE")
#@ String	overlayChoice			(label = "Nuclei/cell overlay choice", choices={"nucleus/cell ID","foci count","both","none"}, description="Select which numbers are added to the overlay.")
#@ ColorRGB	fontColorCells			(label = "Nuclei/cell label color", value="orange", description="The color of the nucleus/cell ID text overlay.")
#@ ColorRGB	fontColorFoci			(label = "Foci count label color", value="red", description="The color of the foci count text overlay.")
#@ Integer	labelFontSize			(label = "Nuclei/cell label font size", value=12, min=4, description="The size of the nuclei/cell ID text overlay.")
#@ String	overlayBrightness		(label = "Nuclei/cell outline brightness", choices={"bright","dim"}, description="The brightness of the nuclei outlines overlay.")
#@ String	outlineColor			(label = "Nuclei/cell outline color", choices={"White","Red","Green","Blue","Cyan","Magenta","Yellow"}, value = "Cyan", description="The color of the nuclei outlines overlay.")
#@ Boolean	debugMode				(label = "Debug mode (show intermediate images)", value=false, description="Used for development and bug fixing: checking this option will trigger displaying many intermediate results during the processing. It will also slow down the analysis.")
#@ String	file_and_image_message0	(value = "<html><p style='font-size:12px'>Need help? Visit the <a href=https://imagej.net/plugins/foci-analyzer>Foci Analyzer</a> website on ImageJ.net</p></html>", visibility="MESSAGE")

version = 1.88;

requires("1.54i");	//Minimum required ImageJ version

/* KNOWN ISSUES
 *  
 * ! Make a possibility to skip the 3D Cellpose parameters dialog
 * ! Foci outside nuclei: Make isotropic - dilate - make non-isotropic
 * ! 3D outlines also in colocalization image
 * 
 * 
 * TO DO | IDEAS
 * 
 * * Open timelapse output files, concatenate, and restore overlays. Now via the script 'concatenate timelapse with overlays'.
 * 
 * * Include 3D Cellpose settings when loading settings from previously analyzed image
 * 
 * * Output:
 * (Centroids of all foci? (MorphoLibJ Analyze regions?))
 * (Distance from edge for all foci?)
 * 
 * * Include Voronoi-Otsu segmentation
 * 
 * * Use Roi.setMinStrokeWidth() for a better scaling of outlines?
 * 
 * * Do not remove lines in log window when areaMaximum detection is used.
 *  
 * * Use ALT/right mouse click to pop up options menu (background subtraction, nuclear segmentation, foci size, etc.
 *   and allow the user to determine the right settings! See MosaicExplorerJ macro.
 *   Currently, pressing ALT during opening the image allows for making a selection.
 * 
 * * Add more: description="..." to the options (onmouseover)
 * 
 * * Background intensity is now very crude (just the mean outside the nuclei). So: mask the background, set (percentile) threshold, measure mean/median. 
 * 
 */

//Check if all dependency Update Sites are installed:
// * - 3D ImageJ Suite
// * - CLIJ
// * - CLIJ2
// * - CLIJx-assistent
// * - CLIJx-assistent-extensions
// * - CSBDeep
// * - IJPB-plugins
// * - PT-BIOP
// * - StarDist


missingPlugin = "";
List.setCommands;
if (List.get("CLIJ2 Macro Extensions")=="") missingPlugin += "CLIJ2, ";
if (List.get("3D Manager")=="") missingPlugin += "3D Image Suite, ";
if (List.get("MorphoLibJ Marker-controlled Watershed (CLIJx, experimental)")=="") missingPlugin += "CLIJx, ";
if (List.get("Intensity Measurements 2D/3D")=="") missingPlugin += "MorphoLibJ (IJPB-Plugins), ";
//if (List.get("Run your network")=="") missingPlugin += "CDBDeep, ";
//if (List.get("Command From Macro")=="") missingPlugin += "StarDist, ";
if (missingPlugin != "") {
	print("\\Clear");
	missingPlugin = missingPlugin.substring(0, missingPlugin.length-2);
	print("Error: Required plugin(s) not found:\n"+missingPlugin+"\n \nFoci Analyzer requires the following Fiji Update Sites to be activated:\n* 3D ImageJ Suite\n* CLIJ\n* CLIJ2\n* CLIJx-assistent\n* CLIJx-assistent-extensions\n* CSBDeep\n* IJPB-plugins\n* StarDist\n \nGo to Help -> Update... -> Manage Update Sites and check the relevant boxes.");
	exit("Error: Required plugin(s) not found:\n"+missingPlugin+"\n \nFoci Analyzer requires the following Fiji Update Sites to be activated:\n* 3D ImageJ Suite\n* CLIJ\n* CLIJ2\n* CLIJx-assistent\n* CLIJx-assistent-extensions\n* CSBDeep\n* IJPB-plugins\n* StarDist\n \nGo to Help -> Update... -> Manage Update Sites and check the relevant boxes.\nThis info is also printed to the Log Window.");
}
List.clear();

if(loadSettingsFromFile == true) {	//TO DO: Add 3D Cellpose settings
	setBatchMode(true);
	loadSettingsPath = File.openDialog("Select a previously analyzed output file (.zip) to load settings from");
	if(endsWith(loadSettingsPath, ".zip")) open(loadSettingsPath);
	else exit("Please select an output .zip file from a previous analysis. Exiting macro.");
	print("Loading settings from "+loadSettingsPath);
	loadVersion = Property.getNumber("00. Version ");
	if(loadVersion < version) showMessage("Warning: The file you have selected has been analyzed with version "+loadVersion+", while you are currently running version "+version+".\nUnexpected things may happen.");
	nucleiChannel = Property.getNumber("02. Nuclei channel ");
	cytoChannel = Property.getNumber("03. Cytoplasm/membrane channel  (-1 if not used) ");
	fociChannelA = Property.getNumber("04. Foci channel A ");
	fociChannelB = Property.getNumber("05. Foci channel B ");
	detectFociChannelB = Property.getNumber("06. Also detect foci channel B? ");
	ThreeDHandling = Property.get("07. 3D handling? ");
	singleSlice = Property.getNumber("08. [single z-slice foci detection only] Slice nr ");
	cropBorder = Property.getNumber("09. Remove image borders (pixels) ");
	XYBinning = Property.getNumber("10. Image XY binning before analysis [1-n] ");
	nucleiSegmentationChoice = Property.get("11. Nuclei/cell segmentation method ");
	downsampleFactorStarDist = Property.getNumber("12. Stardist nuclei rescaling factor [1-n], 0 for automatic, 1 for no rescaling ");
	probabilityThreshold = Property.getNumber("13. Probablility/flow threshold [0.0-1.0] (StarDist/Cellpose) ");
	CellposeModel = Property.getNumber("14. Cellpose model ");
	CellposeDiameter = Property.getNumber("15. Cellpose cell diameter (pixels), 0 for automatic ");
	minNucleusSize_setting = Property.getNumber("16. Remove nulei/cells with diameter smaller than (µm) ");
	maxNucleusSize_setting = Property.getNumber("17. Remove nulei/cells with diameter larger than (µm) ");
	manualRemoveNuclei = Property.get("18. Manually remove segmentated nuclei/cells");
	excludeOnEdges = Property.getNumber("19. Exclude nulei on image edges ");
	//optimizationMode = Property.getNumber("20. Enable foci parameters optimization mode? ");
	fociSizeA = Property.get("21. Foci size channel A ");
	fociSizeB = Property.get("22. Foci size channel B ");
	fociDetectionMethod = Property.get("23. Foci detection method ");
	thresholdFactorA = Property.getNumber("24. Foci intensity threshold bias channel A ");
	thresholdFactorB = Property.getNumber("25. Foci intensity threshold bias channel B ");
	//thresholdA = Property.getNumber("26. Actual threshold value channel A ");
	//thresholdB = Property.getNumber("27. Actual threshold value channel B ");
	minFociSize = Property.getNumber("28. Minimum foci size  (area)");
	maxFociSize = Property.getNumber("29. Maximum foci size  (area)");
	maxFociDistanceOutsideNuclei_setting = Property.getNumber("30. Max distance of foci outside nuclei/cells (µm); -1 for full image");
	minOverlapSize = Property.getNumber("31. Minimum overlap of foci to colocalize (pixels/voxels) ");
	overlayChoice = Property.get("32. Nuclei overlay choice ");
	fontColorCells = Property.get("33. Nuclei label color ");
	fontColorFoci = Property.get("34. Foci count label color ");
	labelFontSize = Property.getNumber("35. Nuclei label size ");
	overlayBrightness = Property.get("36. Nuclei/cell outline bightness ");
	outlineColor = Property.get("37. Nuclei outline color ");
	close();
	setBatchMode(false);
}


//nuclei detection
nucleiBlurRadiusXY = 2;
nucleiBlurRadiusZ = 2;
nucleiMedian3DradiusXY = 2;
nucleiMedian3DradiusZ = 2;
//minNucleusSize = 4;
//maxNucleusSize = 40;
excludeOnEdges = excludeOnEdges;
StarDistDownsampleInterpolation = "Bilinear";	//None, Bilinear or Bicubic
if(nucleiSegmentationChoice == "Cellpose segmentation 2D (or on 2D projection)" || nucleiSegmentationChoice=="Cellpose segmentation 3D") cellpose = true;
else cellpose=false;

//foci detection
var fociFilterSizeXY;	//foci radius for DifferenceOfGaussians
var fociFilterSizeZ;
var fociSizeXY;			//foci radius for detectMaximaBox
var fociSizeZ;

maxTileSize = 2048;			// Maximum StarDist tile size
thresholdMultiplier = 3;	//Default threshold multiplier - foci should be this number of times more intense than the standard deviation of the background (average of all nuclei) 
minFociSize = minFociSize;
maxFociSize = maxFociSize;
var flipSlicesAndChannels = false;
//H_max_peak_flooding = 100;	//Percentage of the peak

thickOutlines = true;		//Width of nuclei outlines, 1 or 2 pixels
thickOutlinesThresholdSize = 1200;	//Above this image size nuclei outlines are always thick (2 pixels)
useLargeSpots = false;		//Make spots in overlay larger
LABELOPACITY = 100;			//opacity of outlines
addNumbersOverlay = true;
saveOverlayImage = true;

var hideCellpose3DDialog = false;
var optimizationMode = optimizationModeSetting;
var processchoice;
var doneOptimizing;
var processAllOtherImages = false;
var processAllOtherImagesFixedThreshold = false;
var useROI = false;
var reAnalyzeFullImage = false;
var currentFileNr = 0;
var displayMode = "composite";	//Default display settings at optimization phase
var activeChannels = "11111";
var maxDisplaySetting;
var processTime = 0;
//var threshold = 0;
var thresholdA = 0;
var thresholdB = 0;
var manualThresholdA = 0;
var manualThresholdB = 0;
var ROIsFolder = "";
var labelImageFolder = "";
var setROIsFolder = false;
var setLabelImageFolder = false;

//Create the azure and orange LUTs for the nuclei and cells (TO DO: maybe do this in a more elegant way than globals)
b_reds = newArray(256);
b_greens = newArray(256); 
b_blues = newArray(256);
create_azure_lut(b_reds, b_greens, b_blues);
o_reds = newArray(256);
o_greens = newArray(256); 
o_blues = newArray(256);
create_orange_lut(o_reds, o_greens, o_blues);

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
if(outputFolder == "0" && useInputAsOutput == false) {
	print("[WARNING] The output folder is not set. Using the input folder as output folder.");
	useInputAsOutput = true;
}
if(!File.exists(outputFolder) && useInputAsOutput == false) {
	createOutputFolder = getBoolean("The output folder '"+outputFolder+"' does not exist. Create?", "Of course, go ahead!", "See if I care!");
	if(createOutputFolder) File.makeDirectory(outputFolder);
	else {
		formatHardDrive = getBoolean("Allright, try this:\nFormat the hard drive?", "Yes, goodbye forever", "No! Please mr. Foci Analyzer, I'll do anything you ask!");
		if(formatHardDrive) {
			showMessage("Ok, you wished for it!");
			exit("Oh, wait. I'll erase myself as well.\nCall it your lucky day then!");
		}
		else exit("That's more like it. Now, run the macro and try again.");
	}
}

if (nucleiSegmentationChoice == "Load ROIs from file" && useInputAsOutput == false && setROIsFolder == false) {
	ROIsFolder = call("ij.Prefs.get", "ROIs.Folder", File.getParent(files[0]));
	Dialog.createNonBlocking("Select a ROIs folder");
	Dialog.addMessage("ROI files should have the same name as the input images without extensions, followed by '_ROIs.zip'.");
	Dialog.addDirectory("Folder containing ROI .zip files", ROIsFolder);
	Dialog.show();
	ROIsFolder = Dialog.getString();
	call("ij.Prefs.set", "ROIs.Folder", ROIsFolder);
	print("Getting segmentations from ROIs in "+ROIsFolder);
	setROIsFolder = true;
}

if (nucleiSegmentationChoice == "Load label images" && useInputAsOutput == false && setLabelImageFolder == false) {
	labelImageFolder = call("ij.Prefs.get", "label.Image.Folder", File.getParent(files[0]));
	Dialog.createNonBlocking("Select a label image folder");
	Dialog.addMessage("Label image files should have the exact same name as the input images.");
	Dialog.addDirectory("Folder containing label image files", labelImageFolder);
	Dialog.show();
	labelImageFolder = Dialog.getString();
	call("ij.Prefs.set", "label.Image.Folder", labelImageFolder);
	print("Getting segmentations from label images in "+labelImageFolder);
	setLabelImageFolder = true;
}
if (cellpose == true && CellposeModel == "custom") {

	CellposeModelPath = call("ij.Prefs.get", "Cellpose.custom.model.path", File.getParent(files[0]));
	CellposeDiameterPref = call("ij.Prefs.get", "Cellpose.diameter", CellposeDiameter);

	Dialog.createNonBlocking("Select a custom Cellpose model");
	Dialog.addFile("Select or drag&drop custom Cellpose model", CellposeModelPath, 80);
	if(CellposeDiameter==0) Dialog.addNumber("Automatic diameter estimation does not work with custom models. Enter a cell diameter", CellposeDiameterPref, 0, 5, "pixels");
	if(cytoChannel>0 && nucleiChannel>0) Dialog.addMessage("N.B. The custom model should be trained on images with a cytoplasm/membrane channel *and* a nucleus channel!\nIf not, no cells will be detected, but Fiji will not crash.");
	else if(cytoChannel>0 && nucleiChannel<=0) Dialog.addMessage("N.B. The custom model should be trained on images with a cytoplasm/membrane channel *only*!\nIf not, no cells will be detected, but Fiji will not crash.");
	else if(cytoChannel<=0 && nucleiChannel>0) Dialog.addMessage("N.B. The custom model should be trained on images with a nucleus channel *only*!\nIf not, no cells will be detected, but Fiji will not crash.");
	Dialog.show();
	CellposeModelPath = Dialog.getString();
	if(CellposeDiameter==0) CellposeDiameter = Dialog.getNumber();
	CellposeModel = CellposeModelPath;		//Both need to be the same

	call("ij.Prefs.set", "Cellpose.custom.model.path", CellposeModelPath);
	call("ij.Prefs.set", "Cellpose.diameter", CellposeDiameter);

	print("Using Cellpose model "+CellposeModelPath);
}
else CellposeModelPath = "path\\to\\own_cellpose_model";	//dummy name


for (currentFileNr = 0; currentFileNr < nrOfImages; currentFileNr++) {
	print("\nProcessing file "+currentFileNr+1+"/"+nrOfImages+": "+files[currentFileNr] + "\n");
	if(processOnlyExtension == "" || endsWith(files[currentFileNr], processOnlyExtension)) processFile(currentFileNr, files[currentFileNr], outputFolder);
	else print("Skipping "+File.getName(files[currentFileNr])+" because it is not a "+processOnlyExtension+" file.");
}
close("Results");

logWindowContents = getInfo("log");
File.saveString(logWindowContents, outputFolder + File.separator + "Log.txt");

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

	if(useInputAsOutput) outputFolder = File.getDirectory(file);
	if(useInputAsOutput) ROIsFolder = File.getDirectory(file);
	if(useInputAsOutput) labelImageFolder = File.getDirectory(file);
	
	
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

	if(nameHasExtension) imageName = File.getNameWithoutExtension(image);
	else imageName = image;
	
	//Initialize image and table
	getDimensions(gwidth, gheight, gchannels, gslices, gframes); // global variables
	getVoxelSize(pixelWidth, pixelHeight, pixelDepth, unit);
	if(bitDepth() == 24) {
		print("[WARNING] The input image is RGB. It will be converted to composite multichannel image with channels 1:Red, 2:Green; 3:Blue.");
		run("Make Composite");
	}
	if(gslices>1) imageIs3D = true;
	else imageIs3D = false;
	
	//timelapse handling - save individual frames and call processFile recursively [smiley with sunglasses]
	if(gframes > 1){
		original = getTitle;
		for(t=1; t<=gframes; t++) {
			showStatus("Splitting time frames... "+t+"/"+gframes);
			showProgress(t, gframes);
			run("Duplicate...", "duplicate frames="+t);
			saveAs("tiff", outputFolder + File.separator + imageName + "__t="+IJ.pad(t, 3));
			close();
		}
		if (nucleiSegmentationChoice == "Load label images") {
			open(labelImageFolder + File.separator + File.getNameWithoutExtension(image) + ".tif");
			for(t=1; t<=gframes; t++) {
				showStatus("Splitting label image time frames... "+t+"/"+gframes);
				showProgress(t, gframes);
				run("Duplicate...", "duplicate range="+t+"-"+t);
				saveAs("tiff", labelImageFolder + File.separator + imageName + "__t="+IJ.pad(t, 3));
				close();
			}
		}
		for(t=1; t<=gframes; t++) {
			processFile(t-1, outputFolder + File.separator + imageName + "__t="+IJ.pad(t, 3)+".tif", outputFolder);
		}
//overlay_image = "Foci_overlay_ch3";
//		for(t=1; t<=gframes; t++) {
//			open(outputFolder + File.separator + imageName + "__t="+IJ.pad(t, 3) + "__" + overlay_image + ".zip");
////			run("To ROI Manager");
////			Overlay.copy;
//			if(t>1) run("Concatenate...", "image1=all_frames image2=" + imageName + "__t="+IJ.pad(t, 3) + "__" + overlay_image + ".tif");
//			rename("all_frames");
//		}
//setBatchMode("exit and display");
		return;
	}

	original = getTitle();
	print("Image path: "+image+"\n");

	if(gchannels == 1 && gslices > 1) {
		if(flipSlicesAndChannels == false) {
			setBatchMode("show");
			Dialog.createNonBlocking("Single channel detected");
			Dialog.addMessage("Warning: ["+imageName+"] has only 1 Channel, but "+gslices+" Slices.");
			Dialog.addRadioButtonGroup("Do you want to flip Channels and Slices?", newArray("Yes, and do this for all subsequent images", "Yes, only for this image", "No (exit macro)"), 3, 1, "Yes, and do this for all subsequent images");
			Dialog.show();
			answer = Dialog.getRadioButton();
			if(answer == "Yes, and do this for all subsequent images") flipSlicesAndChannels = true;
			if(answer == "Yes, only for this image") run("Re-order Hyperstack ...", "channels=[Slices (z)] slices=[Channels (c)] frames=[Frames (t)]");
			else if(answer == "No (exit macro)") exit();
		}
		if(flipSlicesAndChannels == true) run("Re-order Hyperstack ...", "channels=[Slices (z)] slices=[Channels (c)] frames=[Frames (t)]");
		setBatchMode("hide");
	}
	Stack.setDisplayMode("grayscale");
	getDimensions(gwidth, gheight, gchannels, gslices, gframes);
	Stack.setChannel(nucleiChannel);

	if(cropBorder>0) {
		makeRectangle(cropBorder, cropBorder, gwidth-2*cropBorder, gheight-2*cropBorder);
		run("Crop");
		getDimensions(gwidth, gheight, gchannels, gslices, gframes); // global variables
	}
	if(gwidth > thickOutlinesThresholdSize && gheight > thickOutlinesThresholdSize) {
		thickOutlines = true;
	}
	
	if(gslices > 1) anisotropyFactor = pixelDepth / pixelWidth;
	else anisotropyFactor = 0;
	bits = bitDepth();
	
//TO DO: FIND A BETTER WAY FOR THIS:
	if(bits == 8) run("16-bit");	//Convert to 16-bit, because 8-bit restricts the foci labelmap to 255 foci

	Stack.setSlice(gslices/2);
	setBatchMode("show");
	run("Enhance Contrast", "saturated=0.35");

	resultTable = "Foci results per cell";
	if(isOpen(resultTable)) Table.reset(resultTable);
	else {
		Table.create(resultTable);
		Table.setLocationAndSize(0, 0, 1000, 500);
	}
	allFociResultsTable = "All foci statistics";
	if(isOpen(allFociResultsTable)) close(allFociResultsTable);
	Table.create(allFociResultsTable);
	Table.setLocationAndSize(0, 500, 1000, 500);

	if(gslices > 1 && ThreeDHandling == "Detect foci on the Maximum Intensity Projection") {
		selectWindow(original);
		run("Z Project...", "projection=[Max Intensity]");
		setBatchMode("show");
		original = getTitle();
		gslices = 1;
	}
	if(gslices > 1 && ThreeDHandling == "Detect foci on the Standard Deviation Projection") {
		selectWindow(original);
		run("Z Project...", "projection=[Standard Deviation]");
		setBatchMode("show");
		original = getTitle();
		gslices = 1;
	}
	if(gslices > 1 && ThreeDHandling == "Detect foci on the Summed Intensity Projection") {
		selectWindow(original);
		run("Z Project...", "projection=[Sum Slices]");
		setBatchMode("show");
		original = getTitle();
		gslices = 1;
		run("Conversions...", "scale");
	}
	if(gslices > 1 && ThreeDHandling == "Detect foci on a Extended Depth of Focus Projection") {	//TO DO: Doesn't work yet, because multichannel is not compatible with CLIJ2
		selectWindow(original);
		getDimensions(width, height, channels, slices, frames);
		mergeString = "";
		for(c=1; c<=channels; c++) {
			selectWindow(original);
			Stack.setChannel(c);
			Ext.CLIJ2_push(original);
			EDF = "EDF";
			Ext.CLIJ2_extendedDepthOfFocusSobelProjection(original, EDF, 10);
			Ext.CLIJ2_pull(EDF);
			rename("EDF_channel_"+c);
			mergeString += " c"+c+"=EDF_channel_"+c;
		}
		Ext.CLIJ2_release(EDF);
		run("Merge Channels...", mergeString+" create");
		rename(original+"_EDF");
		Stack.setDisplayMode("grayscale");
		original = getTitle();
		setBatchMode("show");
		gslices = 1;
	}
	if(gslices == 1) pixelDepth = 1;
	
	//Run analysis only on a ROI
	wait(50);
	x_ROI = 0;
	y_ROI = 0;
	useROI = false;
	if(isKeyDown("alt")) {
		setTool("rectangle");
		waitForUser("Select a (rectangular) ROI to optimize the analysis on");
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
	if(useROI == false) process_image = original;

	//Segment and label nuclei
	if (nucleiSegmentationChoice == "StarDist nuclei segmentation 2D (or on 2D projection)") nuclei_info = segmentNucleiStarDist(process_image, nucleiChannel, probabilityThreshold, pixelWidth, unit, resultTable);
	else if (nucleiSegmentationChoice == "Classic segmentation") nuclei_info = segmentNucleiClassic(process_image, nucleiChannel, nucleiBlurRadiusXY, nucleiBlurRadiusZ, nucleiMedian3DradiusXY, nucleiMedian3DradiusZ);
	else if (cellpose == true) 									 nuclei_info = segmentCellsCellpose(process_image, nucleiChannel, cytoChannel, probabilityThreshold, pixelWidth, pixelDepth, unit, resultTable);
	else if (nucleiSegmentationChoice == "Load ROIs from file")  nuclei_info = loadROIs(original, nucleiChannel, ROIsFolder);
	else if (nucleiSegmentationChoice == "Load label images") 	 nuclei_info = loadROIs(original, nucleiChannel, labelImageFolder);

	if(nuclei_info[0] == "FileNotFound") continue;
	labelmap_nuclei = nuclei_info[0];
	nuclei_outlines = nuclei_info[1];
	nrNuclei = nuclei_info[2];

	if(manualRemoveNuclei != "No thanks" && nrNuclei > 0) {
		nuclei_info = manually_remove_labels(labelmap_nuclei, nuclei_outlines, nrNuclei, process_image, imageName);
		labelmap_nuclei = nuclei_info[0];
		nuclei_outlines = nuclei_info[1];
		nrNuclei = nuclei_info[2];
//		labelmap_nuclei = "Labelmap_nuclei_edited";
	}
	if(nrNuclei > 0) {
		//Write ID and area to foci results table
		run("Clear Results");
		Ext.CLIJ2_statisticsOfLabelledPixels(labelmap_nuclei, labelmap_nuclei);
		nucleus_id_ = Table.getColumn("IDENTIFIER", "Results");
		nucleus_area_ = Table.getColumn("PIXEL_COUNT", "Results");
		nucleus_area_ = multiplyArraywithScalar(nucleus_area_, Math.sqr(pixelWidth));
		Table.setColumn("Cell ID", nucleus_id_, resultTable);
		if (nucleiSegmentationChoice == "Load ROIs from file") {
			for (i = 0; i < nrNuclei; i++) {
				roiManager("select", i);
				Table.set("Cell UUID", i, Roi.getName, resultTable);
			}
			roiManager("deselect");
		}
		Table.setColumn("Cell area 2D ("+unit+"^2)", nucleus_area_, resultTable);

		//Create a 3D version of the 2D nuclei labelmap, if required
		Ext.CLIJ2_getDimensions(labelmap_nuclei, labelmap_width, labelmap_height, labelmap_depth);
		if(gslices>1 && labelmap_depth==1) Ext.CLIJ2_imageToStack(labelmap_nuclei, labelmap_nuclei_3D, gslices);
		else labelmap_nuclei_3D = labelmap_nuclei;

		//Foci filtering and detection - in a loop to enable parameter optimization 
		firstTimeProcessing = true;
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
			if(processAllOtherImages == false) doneOptimizing = false;	//Reset this parameter from the previous round
			
			//Create dilated labelmap for including foci outside the nuclei
			labelmap_nuclei_3D_dilated = "labelmap_nuclei_3D_dilated";
			if(maxFociDistanceOutsideNuclei_setting != 0) {
				maxFociDistanceOutsideNuclei = round(maxFociDistanceOutsideNuclei_setting / pixelWidth);
				if(maxFociDistanceOutsideNuclei_setting < 0) Ext.CLIJ2_dilateLabels(labelmap_nuclei_3D, labelmap_nuclei_3D_dilated, maxOf(gwidth, gheight));
				else Ext.CLIJ2_dilateLabels(labelmap_nuclei_3D, labelmap_nuclei_3D_dilated, maxFociDistanceOutsideNuclei);
				if(ThreeDHandling == "Detect foci in 3D (or 2D when N/A)" && nucleiSegmentationChoice == "Cellpose segmentation 3D" && pixelDepth/pixelWidth > 1.33) print("[WARNING] Due to non-isotropic pixels the foci region will be expanded from the nuclei more in Z than in X and Y! (by a factor of "+d2s(pixelDepth/pixelWidth,1)+")");
				//N.B. For non-isotropic 3D data this 3D dilation is not fair, but making the labelmap isotropic creates intermediate (non-integer) values. Oh well..
			}
			else Ext.CLIJ2_copy(labelmap_nuclei_3D, labelmap_nuclei_3D_dilated);
			if(debugMode) showImagefromGPU(labelmap_nuclei_3D_dilated);
	
			//Create outlines from dilated nuclei labelmap
			nuclei_dilated_outlines = "nuclei_dilated_outlines";
			if(isOpen("nuclei_dilated_outlines")) close(nuclei_dilated_outlines);
			if(maxFociDistanceOutsideNuclei_setting != 0 && nucleiSegmentationChoice != "Cellpose segmentation 3D") {
				Ext.CLIJ2_copySlice(labelmap_nuclei_3D_dilated, labelmap_nuclei_2D_dilated, 0);
				Ext.CLIJ2_detectLabelEdges(labelmap_nuclei_2D_dilated, nuclei_dilated_edges);
				Ext.CLIJ2_mask(labelmap_nuclei_2D_dilated, nuclei_dilated_edges, nuclei_dilated_outlines);
				Ext.CLIJ2_release(nuclei_dilated_edges);
				Ext.CLIJ2_release(labelmap_nuclei_2D_dilated);
				Ext.CLIJ2_pullBinary(nuclei_dilated_outlines);
				if(outlineColor == "White") outlineColor = "Grays";
				run(outlineColor);
			}
			else if(maxFociDistanceOutsideNuclei_setting != 0 && nucleiSegmentationChoice == "Cellpose segmentation 3D") {
				Ext.CLIJ2_maximumSliceBySliceSphere(labelmap_nuclei_3D_dilated, labelmap_maximum, 1, 1);
				Ext.CLIJ2_minimumSliceBySliceSphere(labelmap_nuclei_3D_dilated, labelmap_minimum, 1, 1);
				Ext.CLIJ2_subtractImages(labelmap_maximum, labelmap_minimum, label_edges);
				Ext.CLIJ2_release(labelmap_maximum);
				Ext.CLIJ2_release(labelmap_minimum);
				Ext.CLIJ2_threshold(label_edges, nuclei_dilated_outlines, 1);
				Ext.CLIJ2_pull(nuclei_dilated_outlines);
				Ext.CLIJ2_release(label_edges);
				Ext.CLIJ2_release(nuclei_dilated_outlines);
				if(outlineColor == "White") outlineColor = "Grays";
				run(outlineColor);
			}

			//Check if fociChannelB actually exists and is different from Foci channel A
			if(detectFociChannelB == true && fociChannelB > gchannels) exit("Error: The selected Foci channel B ("+fociChannelB+") is higher than the number of channels of the image ("+gchannels+")."); 
			if(detectFociChannelB == true && fociChannelB == fociChannelA) {
				showMessageWithCancel("Selected Foci channels are the same", "Warning: The selected Foci channel B ("+fociChannelB+") is the same as Foci channel A!\nIf you continue only this channel will be analyzed.");
				detectFociChannelB = false;
			}

			//Foci channel A
			detections_fociA = detect_foci(process_image, fociChannelA, fociSizeA, anisotropyFactor, firstTimeProcessing, labelmap_nuclei, labelmap_nuclei_3D, thresholdFactorA, thresholdA, "A");
			labelmap_fociA = detections_fociA[0];
			mask_fociA = detections_fociA[1];
			spots_fociA = detections_fociA[2];
			thresholdA = detections_fociA[3];

			if(firstTimeProcessing == false) call("ij.gui.ImageWindow.setNextLocation", x_image, y_image);
			overlayA = mergeOriginalAndDetection(original, nrNuclei, nuclei_outlines, mask_fociA, spots_fociA, fociChannelA, useROI, x_ROI, y_ROI);
			overlay_image = overlayA;	//Will be overwritten if channel B is also used
		
			if(firstTimeProcessing == false && detectFociChannelB == false) close("Processing...");
			//Foci channel B
			if(detectFociChannelB == true) {
				detections_fociB = detect_foci(process_image, fociChannelB, fociSizeB, anisotropyFactor, firstTimeProcessing, labelmap_nuclei, labelmap_nuclei_3D, thresholdFactorB, thresholdB, "B");
				labelmap_fociB = detections_fociB[0];
				mask_fociB = detections_fociB[1];
				spots_fociB = detections_fociB[2];
				thresholdB = detections_fociB[3];

				overlayB = mergeOriginalAndDetection(original, nrNuclei, nuclei_outlines, mask_fociB, spots_fociB, fociChannelB, useROI, x_ROI, y_ROI);
				if(firstTimeProcessing == false) close("Processing...");
				Overlay.copy();	//Preserve nuclei outline overlays
				overlay_image = "Foci_overlay_ch"+fociChannelA+"_and_ch"+fociChannelB;
				run("Concatenate...", "  title="+overlay_image+" image1="+overlayA+" image2="+overlayB+" image3=[-- None --]");
				run("Stack to Hyperstack...", "order=xyczt(default) channels="+nSlices/(gslices*2)+" slices="+gslices+" frames=2 display=Composite");
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
				if(cytoChannel>0) {
					Stack.setChannel(5);
					setMinAndMax(minDisplayCells, maxDisplayCells);
				}
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
			selectWindow("Channels");
			setLocation(0, 0);
			run("Brightness/Contrast...");
			selectWindow("B&C");
			setLocation(-5, 210);
			if(firstTimeProcessing == true) getLocationAndSize(x_image, y_image, imageWidth, imageHeight);
			//if(slices>1) run("Animation Options...", "speed=4 loop start");
			
			//Parameter optimization dialog
			if(optimizationMode == true && processAllOtherImages == false) {

				Dialog.createNonBlocking("Optimize settings for foci detection");
				Dialog.addMessage("Inspect the detected foci and optimize the settings. Some tips:\n• Hide and show the foci with the checkbox 'channel 3' in the Channels Tool (upper left corner of the screen).\n• Change the display brighness of the 'channels' (1:nuclei, 2:foci, 3:detected foci, 4:foci centers, [5:cells]) in the B&C window.\n• Zoom [+/- or up/down keys] and pan [hold space & drag the image].\n• Use the sliders below the image to change the active channel [arrow keys] and displayed z-slice [Ctrl + arrow keys],\n   and to switch between Foci channels A and B [Alt + arrow keys].\n\n ", 12, "#000080");
				Dialog.addChoice("Foci size channel "+fociChannelA, newArray("tiny (0.5 pixels)","small (1 pixel)","average (2 pixels)","large (4 pixels)","huge (8 pixels)","other (define later)"), fociSizeA);
				if(detectFociChannelB) Dialog.addChoice("Foci size channel "+fociChannelB, newArray("tiny (0.5 pixels)","small (1 pixel)","average (2 pixels)","large (4 pixels)","huge (8 pixels)","other (define later)"), fociSizeB);
				Dialog.addChoice("Detection method", newArray("Marker-controlled watershed (recommended)","AreaMaxima local maximum detection", "Manual thresholding"), fociDetectionMethod);
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
				Dialog.addNumber("Extend foci detection outside nuclei/cells with", maxFociDistanceOutsideNuclei_setting, 1, 4, unit+" (-1 for full image)")
				if(detectFociChannelB) {
					if(gslices > 1) Dialog.addNumber("Minimum overlap of foci", minOverlapSize, 1, 4, "voxels");
					else Dialog.addNumber("Minimum overlap of foci", minOverlapSize, 1, 4, "pixels");
				}
				Dialog.addChoice("Nuclei outline visualization", newArray("bright","dim"), overlayBrightness);
				items = newArray("Recalculate with these settings", "Done optimizing | Process and optimize the next image", "Done optimizing | Process all images with these settings (calculate thresholds for each image)", "Done optimizing | Process all images with these settings (fix current threshold levels)");
				Dialog.addChoice("Action", items, "Recalculate with these settings");
				Dialog.addMessage("Click Help for more info (Foci Analyzer ImageJ site)", 12, "#000080");
				Dialog.addHelp("https://imagej.net/plugins/foci-analyzer");

				//determine dialog location
				if(x_image + imageWidth + 500 < screenWidth) Dialog.setLocation(x_image+imageWidth-15, y_image);
				else if(x_image > 500) Dialog.setLocation(x_image-530, y_image);
//				else if(y_image + gheight + 300 < screenHeight) Dialog.setLocation(x_image, y_image+gheight+100);
				else Dialog.setLocation(x_image, y_image);
				Dialog.show();

				//Get dialog entries
				fociSizeA = Dialog.getChoice();
				if(detectFociChannelB == true) fociSizeB = Dialog.getChoice();
				fociDetectionMethod = Dialog.getChoice();
				thresholdFactorA = Dialog.getNumber();
				if(detectFociChannelB) thresholdFactorB = Dialog.getNumber();
				minFociSize = Dialog.getNumber();
				maxFociSize = Dialog.getNumber();
				maxFociDistanceOutsideNuclei_setting = Dialog.getNumber();
				if(detectFociChannelB) minOverlapSize = Dialog.getNumber();
				overlayBrightness = Dialog.getChoice();
				processChoice = Dialog.getChoice();
				if(processChoice == "Done optimizing | Process and optimize the next image") doneOptimizing = true; processAllOtherImages = false;
				if(processChoice == "Done optimizing | Process all images with these settings (calculate thresholds for each image)") { doneOptimizing = true; processAllOtherImages = true; }
				if(processChoice == "Done optimizing | Process all images with these settings (fix current threshold levels)") { doneOptimizing = true; processAllOtherImages = true; processAllOtherImagesFixedThreshold = true; }

				//Write persistent Script parameters for the next time the macro is run.
				setScriptParameterValue("fociSizeA", fociSizeA);
				setScriptParameterValue("fociSizeB", fociSizeB);
				setScriptParameterValue("thresholdFactorA", thresholdFactorA);
				setScriptParameterValue("thresholdFactorB", thresholdFactorB);
				setScriptParameterValue("fociDetectionMethod",fociDetectionMethod);
				setScriptParameterValue("minFociSize",minFociSize);
				setScriptParameterValue("maxFociSize",maxFociSize);
				setScriptParameterValue("maxFociDistanceOutsideNuclei_setting",maxFociDistanceOutsideNuclei_setting);
				setScriptParameterValue("minOverlapSize",minOverlapSize);
				setScriptParameterValue("overlayBrightness",overlayBrightness);
				
				//get image display properties
				selectWindow(overlay_image);
				getLocationAndSize(x_image, y_image, imageWidth, imageHeight);
				Stack.getActiveChannels(activeChannels);
				Stack.getDisplayMode(displayMode);
				Stack.getPosition(currentChannel, currentSlice, currentFrame);
				Stack.setChannel(1);
				getMinAndMax(minDisplayNuclei, maxDisplayNuclei);
				Stack.setChannel(2);
				getMinAndMax(minDisplayFoci, maxDisplayFoci);
				if(cytoChannel>0) {
					Stack.setChannel(5);
					getMinAndMax(minDisplayCells, maxDisplayCells);
				}
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
				else useROI = false;
// N.B. REMOVING THE PREVIOUS LINE WILL RESULT IN THE FULL IMAGE BEING OPENED AND ANALYZED! THIS HAPPENS ANYWAY WHEN OPTIMIZATION IS OFF.

			}
		firstTimeProcessing = false;
		} while(optimizationMode == true && doneOptimizing == false);

		if(useROI == false) {
			
			//Get nuclei positions
			if(maxFociDistanceOutsideNuclei_setting != 0) {
				run("Clear Results");
				Ext.CLIJ2_statisticsOfLabelledPixels(labelmap_nuclei_3D_dilated, labelmap_nuclei_3D_dilated);
			}
			selectWindow("Results");
			boundingBox_X = Table.getColumn("BOUNDING_BOX_X");
			boundingBox_Y = Table.getColumn("BOUNDING_BOX_Y");
			boundingBox_Z = Table.getColumn("BOUNDING_BOX_Z");
			boundingBox_width = Table.getColumn("BOUNDING_BOX_WIDTH");
			boundingBox_height = Table.getColumn("BOUNDING_BOX_HEIGHT");
			boundingBox_depth = Table.getColumn("BOUNDING_BOX_DEPTH");
			Array.getStatistics(boundingBox_width, minWidth, maxWidth);
			Array.getStatistics(boundingBox_height, minHeight, maxHeight);
			Array.getStatistics(boundingBox_depth, minDepth, maxDepth);
			if(debugMode) print("\nMaximum nucleus bounding box: "+maxWidth+", "+maxHeight);		

			//Measure the foci 
			nrFoci = measureFoci(original, fociChannelA, nrNuclei, labelmap_nuclei_3D, labelmap_fociA, boundingBox_X, boundingBox_Y, boundingBox_Z, maxWidth, maxHeight, maxDepth);
			if(detectFociChannelB) nrFoci = measureFoci(original, fociChannelB, nrNuclei, labelmap_nuclei_3D, labelmap_fociB, boundingBox_X, boundingBox_Y, boundingBox_Z, maxWidth, maxHeight, maxDepth);

			if(gslices > 1 && nrFoci>0) Table.renameColumn("Volume", "Volume (voxels)", allFociResultsTable);
			else if (gslices == 1 && nrFoci>0) Table.renameColumn("Volume", "Area (pixels)", allFociResultsTable);
			Table.update(allFociResultsTable);

			Table.deleteRows(0, 0, "Results");	//remove background label
			overlay_numbers_on_image(overlay_image);

	//		print("\n\nGPU Memory after channel "+fociChannelB);
	//		Ext.CLIJ2_reportMemory();
		
			//Compute A->B foci colocalization
			if(detectFociChannelB) {
				foci_overlap_map = computeOverlap("Mask_foci_ch" + fociChannelA, "Mask_foci_ch" + fociChannelB, nrNuclei, labelmap_nuclei_3D, boundingBox_X, boundingBox_Y, maxWidth, maxHeight);	//These images are still open in RAM/GPU
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
			List.set("00. ==Foci Analyzer settings=", "");
			List.set("00. Date ", " " + DayNames[dayOfWeek] + " " + dayOfMonth + " " + MonthNames[month] + " " + year);
			List.set("00. Time ", " " + hour +":"+IJ.pad(minute,2)+":"+IJ.pad(second,2));
			List.set("00. Version ", version);
			if(loadSettingsFromFile == true) List.set("01. Analysis settings were loaded from ", loadSettingsPath);
			else List.set("01. Settins were entered from dialog", 1);
			List.set("02. Nuclei channel ", nucleiChannel);
			List.set("03. Cytoplasm/membrane channel  (-1 if not used) ", cytoChannel);
			List.set("04. Foci channel A ", fociChannelA);
			List.set("05. Foci channel B ", fociChannelB);
			List.set("06. Also detect foci channel B? ", detectFociChannelB);
			List.set("07. 3D handling? ", ThreeDHandling);
			List.set("08. [single z-slice foci detection only] Slice nr ", singleSlice);
			List.set("09. Remove image borders (pixels) ", cropBorder);
			List.set("10. Image XY binning before analysis [1-n] ", XYBinning);
			List.set("11. Nuclei/cell segmentation method ", nucleiSegmentationChoice);
			List.set("12. Stardist nuclei rescaling factor [1-n], 0 for automatic, 1 for no rescaling ", downsampleFactorStarDist);
			List.set("13. Probablility/flow threshold [0.0-1.0] (StarDist/Cellpose) ", probabilityThreshold);
			List.set("14. Cellpose model ", CellposeModel);
			List.set("15. Cellpose cell diameter (pixels), 0 for automatic ", CellposeDiameter);
			List.set("16. Remove nulei/cells with diameter smaller than (µm) ", minNucleusSize_setting);
			List.set("17. Remove nulei/cells with diameter larger than (µm) ", maxNucleusSize_setting);
			List.set("18. Manually remove segmentated nuclei/cells", manualRemoveNuclei);
			List.set("19. Exclude nulei on image edges ", excludeOnEdges);
			List.set("20. Enable foci parameters optimization mode? ", optimizationMode);
			List.set("21. Foci size channel A ", fociSizeA);
			List.set("22. Foci size channel B ", fociSizeB);
			List.set("23. Foci detection method ", fociDetectionMethod);
			List.set("24. Foci intensity threshold bias channel A ", thresholdFactorA);
			List.set("25. Foci intensity threshold bias channel B ", thresholdFactorB);
			List.set("26. Actual threshold value channel A ", thresholdA);
			List.set("27. Actual threshold value channel B ", thresholdB);
			List.set("28. Minimum foci size (area) ", minFociSize);
			List.set("29. Maximum foci size (area) ", maxFociSize);
			List.set("30. Max distance of foci outside nuclei/cells (µm); -1 for full image ", maxFociDistanceOutsideNuclei_setting);
			List.set("31. Minimum overlap of foci to colocalize (pixels/voxels) ", minOverlapSize);
			List.set("32. Nuclei overlay choice ", overlayChoice);
			List.set("33. Nuclei/cell label color ", fontColorCells);
			List.set("34. Nuclei/cell label color ", fontColorFoci);
			List.set("35. Nuclei/cell label size ", labelFontSize);
			List.set("36. Nuclei/cell outline bightness ", overlayBrightness);
			List.set("37. Nuclei/cell outline color ", outlineColor);

			list = List.getList();

			selectWindow(resultTable);
			Table.save(outputFolder + File.separator + imageName + "__Foci_results.tsv");
			selectWindow(allFociResultsTable);
			Table.save(outputFolder + File.separator + imageName + "__All_Foci_statistics.tsv");

			if(saveOverlayImage == true) {
				if(detectFociChannelB == true) {
					selectImage(foci_overlap_map);
					Stack.setSlice(gslices/2);
					setBatchMode("show");
					setMetadata("info", list);
					saveAs("zip", outputFolder + File.separator + imageName + "__foci_coloc_map");
				}
			}
			selectImage(overlay_image);
			setBatchMode("hide");
			setBatchMode("show");	//Lame trick to move this window to the front 
			setMetadata("info", list);
			saveAs("zip", outputFolder + File.separator + imageName + "__" + overlay_image);
			overlay_image = getTitle(); //Name has changed
			run("Clear Results");

//			Ext.CLIJ2_release(labelmap_nuclei);
//			Ext.CLIJ2_release(labelmap_nuclei_3D);
//			Ext.CLIJ2_release(labelmap_nuclei_3D_dilated);
//			Ext.CLIJ2_release("foci_ch"+fociChannelA);
//			if(detectFociChannelB == true) Ext.CLIJ2_release("foci_ch"+fociChannelB);
//			Ext.CLIJ2_release("Labelmap_detected_foci_ch"+fociChannelA);
//			if(detectFociChannelB == true) Ext.CLIJ2_release("Labelmap_detected_foci_ch"+fociChannelB);
//			Ext.CLIJ2_reportMemory();
			Ext.CLIJ2_clear();
		}
	}
	else print("[WARNING] 0 nuclei detected in "+ image +" !");
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

	if(unit == "µm" || unit == "um" || unit == "microns" || unit == "micron") nucleiBlurRadiusXY = 0.5/pixelWidth;	//Set sigma to 0.5 µm
	Ext.CLIJ2_gaussianBlur3D(nuclei, nuclei_filtered1, nucleiBlurRadiusXY, nucleiBlurRadiusXY, nucleiBlurRadiusZ);
	Ext.CLIJ2_maximumZProjection(nuclei_filtered1, nuclei_filtered_MAX);
	Ext.CLIJ2_automaticThreshold(nuclei_filtered_MAX, thresholded, "Otsu");
	Ext.CLIJ2_pullBinary(thresholded);
	run("Watershed");	//Not on GPU because results are not so good.
	if(debugMode) { setBatchMode("show"); roiManager("Show None"); }
	run("Properties...", "unit=&unit pixel_width=&pixelWidth pixel_height=&pixelHeight voxel_depth=1.0000");
	if(excludeOnEdges) run("Analyze Particles...", "size=" + minNucleusSize + "-" + maxNucleusSize + " pixel circularity=0.25-1.00 show=[Count Masks] exclude include add");
	else 			   run("Analyze Particles...", "size=" + minNucleusSize + "-" + maxNucleusSize + " pixel circularity=0.25-1.00 show=[Count Masks] include add");

	rename("Labelmap_nuclei");
	if(debugMode) { run("glasbey on dark"); setBatchMode("show"); }
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
	//Check whether the StarDist dependencies have been installed
	List.setCommands;
	if (List.get("Run your network")=="") missingPlugin += "CDBDeep, ";
	if (List.get("Command From Macro")=="") missingPlugin += "StarDist, ";
	if (missingPlugin != "") {
		print("\\Clear");
		missingPlugin = missingPlugin.substring(0, missingPlugin.length-2);
		print("Error: Required plugin(s) not found:\n"+missingPlugin+"\n \nFoci Analyzer requires the following Fiji Update Sites to be activated:\n* 3D ImageJ Suite\n* CLIJ\n* CLIJ2\n* CLIJx-assistent\n* CLIJx-assistent-extensions\n* CSBDeep\n* IJPB-plugins\n* StarDist\n \nGo to Help -> Update... -> Manage Update Sites and check the relevant boxes.");
		exit("Error: Required plugin(s) not found:\n"+missingPlugin+"\n \nFoci Analyzer requires the following Fiji Update Sites to be activated:\n* 3D ImageJ Suite\n* CLIJ\n* CLIJ2\n* CLIJx-assistent\n* CLIJx-assistent-extensions\n* CSBDeep\n* IJPB-plugins\n* StarDist\n \nGo to Help -> Update... -> Manage Update Sites and check the relevant boxes.\nThis info is also printed to the Log Window.");
	}
	List.clear();

	selectWindow(image);
	getDimensions(width, height, channels, slices, frames);
	minNucleusSize = PI*Math.sqr((minNucleusSize_setting / pixelWidth / 2));	//Calculate the nucleus area as if it were a circle
	maxNucleusSize = PI*Math.sqr((maxNucleusSize_setting / pixelWidth / 2));

	run("Duplicate...", "title=nuclei duplicate channels=" +channel);	//Get nucleus channel
	starDist_input_image = getTitle();

	//Downsample and Z-project.**** TO DO (maybe): Perform this on the GPU (but uses more GPU RAM) and requires splitting and combining the channels.
	if(downsampleFactorStarDist == 0 && unit == "µm" || unit == "um" || unit == "microns" || unit == "micron") {
		downsampleFactorStarDist = 0.25/pixelWidth;	//scale to 0.25 um/pixel
		if(downsampleFactorStarDist < 1 &&  downsampleFactorStarDist > 0.5) downsampleFactorStarDist = 1;	//Do not upsample unless the pixel size is > 0.5 um.
		else if(downsampleFactorStarDist < 0.25 && XYBinning>1) print ("[WARNING] The pixel size is very large ("+pixelWidth+" "+unit+"). Stardist may have difficulties segmenting the nuclei. Try setting XY binning to a lower number (currently "+XYBinning);
		else if(downsampleFactorStarDist < 0.25 && XYBinning==1) print ("[WARNING] The pixel size is very large ("+pixelWidth+" "+unit+"). Stardist may have difficulties segmenting the nuclei.");
		else if(downsampleFactorStarDist > 0.8 && downsampleFactorStarDist < 1.2) downsampleFactorStarDist = 1;	//Too small change - skip rescaling 
		print("Pixel size (after "+XYBinning+"x"+XYBinning+" XY binning): "+pixelWidth+" µm\nStarDist downsample factor (automatic): "+downsampleFactorStarDist+" (effective pixel size: "+pixelWidth*downsampleFactorStarDist+" µm)\n");
	}
	else if (downsampleFactorStarDist == 0) {
		print("[WARNING] Pixel size seems incorrect ("+pixelWidth+" "+unit+")"+". Cannot determine the nuclei downsample factor for Stardist. Stardist may have difficulties segmenting the nuclei.");
		downsampleFactorStarDist = 1;
	}
	else {
		if(downsampleFactorStarDist < 0.15) {
			print("[WARNING] Overruling StarDist downsample factor ("+downsampleFactorStarDist+"). Will be set to the minimum value of 0.15");
			downsampleFactorStarDist = 0.15;	//Allow upscaling up to 6.667x
		}
		print("Pixel size (after "+XYBinning+"x"+XYBinning+" XY binning): "+pixelWidth+" µm\nStarDist downsample factor set to: "+downsampleFactorStarDist+" (effective pixel size: "+pixelWidth*downsampleFactorStarDist+" µm)\n");
	}

	if(downsampleFactorStarDist != 1) {
		run("Scale...", "x="+1/downsampleFactorStarDist+" y="+1/downsampleFactorStarDist+" interpolation="+StarDistDownsampleInterpolation+" average process create");
		rename("nuclei_downscaled");

		starDist_input_image = getTitle();
	}
	if(slices>1) {
		run("Z Project...", "projection=[Max Intensity]");
		rename("nuclei_Zprojected");
		starDist_input_image = getTitle();
		if(downsampleFactorStarDist != 1) close("nuclei_downscaled");
	}

	//Run StarDist
	getDimensions(dswidth, dsheight, dschannels, dsslices, dsframes);
	starDistTiles = pow(floor((maxOf(dswidth, dsheight)/maxTileSize)-1)+1,2);	//Determine the nr. of tiles
	run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':'"+starDist_input_image+"', 'modelChoice':'Versatile (fluorescent nuclei)', 'normalizeInput':'true', 'percentileBottom':'1.0', 'percentileTop':'99.60000000000001', 'probThresh':'"+probabilityThreshold+"', 'nmsThresh':'0.3', 'outputType':'ROI Manager', 'nTiles':'"+starDistTiles+"', 'excludeBoundary':'2', 'roiPosition':'Stack', 'verbose':'false', 'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");

	//Scale up ROIs
	if(downsampleFactorStarDist != 1) RoiManager.scale(downsampleFactorStarDist, downsampleFactorStarDist, false);

	//Spline fit ROIs
	selectWindow(image);	//Have to select an image with the correct size
	for (i = 0; i < roiManager("count"); i++) {
		roiManager("select", i);
		run("Fit Spline");
//		Roi.getSplineAnchors(x, y);
//		Roi.setPolygonSplineAnchors(x, y);
		roiManager("Update");
	}
	
	//Convert ROIs to label map
	run("Select None");
	labelmap_nuclei = "Labelmap_nuclei_unfiltered";
//	run("ROI Manager to LabelMap(2D)");
	labelmap_nuclei = ROI_Manager_to_labelmap(image);
	
	setBatchMode("show");
	run("glasbey_on_dark");
	setMinAndMax(0, 255);

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
	Ext.CLIJ2_getMaximumOfAllPixels(labelmap_nuclei_edges_excluded, nrNucleiBeforeSizeFiltering);
	Ext.CLIJ2_excludeLabelsOutsideSizeRange(labelmap_nuclei_edges_excluded, labelmap_nuclei_filtered, minNucleusSize, maxNucleusSize);
	Ext.CLIJ2_release(labelmap_nuclei_edges_excluded);
	labelmap_nuclei_final = "Labelmap_nuclei";
	Ext.CLIJ2_closeIndexGapsInLabelMap(labelmap_nuclei_filtered, labelmap_nuclei_gapsClosed);
//	if(downsampleFactorStarDist > 1) {
		Ext.CLIJ2_greyscaleOpeningSphere(labelmap_nuclei_gapsClosed, labelmap_nuclei_final, round(downsampleFactorStarDist + 1), round(downsampleFactorStarDist + 1), 0);	//Smooth labels a bit
		Ext.CLIJ2_release(labelmap_nuclei_gapsClosed);
//	}
//	else labelmap_nuclei_final = labelmap_nuclei_gapsClosed;
	Ext.CLIJ2_release(labelmap_nuclei_filtered);
	Ext.CLIJ2_getMaximumOfAllPixels(labelmap_nuclei_final, nrNuclei);	//Count the number of nuclei
	print(nrNucleiBeforeSizeFiltering - nrNuclei + " nuclei were removed due to size restrictions (" + minNucleusSize_setting + " - " + maxNucleusSize_setting + " " + unit + " | "+d2s(minNucleusSize,0)+" - "+d2s(maxNucleusSize,0)+" pixels).");

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


function loadROIs(image, channel, folder) {
//	File.setDefaultDir(outputFolder);
	run("Duplicate...", "title=nuclei duplicate channels=" +channel);	//Get nucleus channel

	roiManager("reset");
	if(nucleiSegmentationChoice == "Load ROIs from file") {
		if(File.exists(folder + File.separator + File.getNameWithoutExtension(image) + "_ROIs.zip")) roiManager("open", folder + File.separator + File.getNameWithoutExtension(image) + "_ROIs.zip");
		else {
			print("WARNING: ROI file not found!  "+folder + File.separator + File.getNameWithoutExtension(image) + "_ROIs.zip\nSkipping this image");
			return newArray("FileNotFound");
		}
		run("ROI Manager to LabelMap(2D)");
	}
	else if(nucleiSegmentationChoice == "Load label images") {
		if(File.exists(folder + File.separator + File.getNameWithoutExtension(image) + ".tif")) open(folder + File.separator + File.getNameWithoutExtension(image) + ".tif");	//Label images must have the same name as the original image
		else {
			print("WARNING: Labelmap file not found!  "+folder + File.separator + File.getNameWithoutExtension(image) + ".tif\nSkipping this image");
			return newArray("FileNotFound");			
		}
	}
	labelmap_nuclei_final = "Labelmap_nuclei";
	rename(labelmap_nuclei_final);
	
	Ext.CLIJ2_push(labelmap_nuclei_final);
//	close(labelmap_nuclei);
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


function segmentCellsCellpose (image, nucleiChannel, cytoChannel, probabilityThreshold, pixelWidth, pixelDepth, unit, resultTable) {
	//Check which version of the BIOP Cellpose wrapper is present, if at all.
	List.setCommands;
	if (List.get("Cellpose ...")!="") {
		//Get Cellpose settings
		//This works for the new wrapper
		envPath = getPref("Packages.ch.epfl.biop.wrappers.cellpose.ij2commands.Cellpose", "env_path");
		envType = getPref("Packages.ch.epfl.biop.wrappers.cellpose.ij2commands.Cellpose", "env_type");
		if(envType == "<null>") envType = "conda";	//Default is conda, but returns <null>
		oldCellposeWrapper = false;
		print("Cellpose environment type: "+envType);
		print("Cellpose environment path: "+envPath);
	}
	else if (List.get("Cellpose Advanced")!="") {
		oldCellposeWrapper = true;
//		This works for the old wrapper, but not for the new one...
//		prefs = split(File.openAsString(getDir("preferences") + "IJ_Prefs.txt"),"\n");
//		for (i = 0; i < prefs.length; i++) {
//			data = split(prefs[i],"=");
//			if(data.length>1) List.set(data[0], data[1]);
//		}
//		envPath = List.get(".ch.epfl.biop.wrappers.cellpose.Cellpose.envDirPath");
//		envPath = envPath.replace("\\\\:", ":");	//This line is needed because sometimes the path has an extra '\' after the drive letter.
//		envType = List.get(".ch.epfl.biop.wrappers.cellpose.Cellpose.envType");
		print("[WARNING] Old version of the Cellpose wrapper detected. It may work now, but may be broken in the future. Please Update Fiji via Help -> Update..."
	}
	else exit("ERROR: Cellpose wrapper not found! Update Fiji and activate the PT-BIOP Update site.");
	List.clear();
	if(nucleiSegmentationChoice == "Cellpose segmentation 3D" && ThreeDHandling != "Detect foci in 3D (or 2D when N/A)" && imageIs3D == true) {
		print("[WARNING] 3D Cellpose segmentation is selected, but '3D image handling' is set to '"+ThreeDHandling+"'.\nProceeding with 2D Cellpose segmentation.");
		nucleiSegmentationChoice = "Cellpose segmentation 2D (or on 2D projection)";
	}
	if(nucleiSegmentationChoice == "Cellpose segmentation 3D" && ThreeDHandling == "Detect foci in 3D (or 2D when N/A)" && imageIs3D == false) {
		print("[WARNING] 3D Cellpose segmentation is selected, but the image has only 1 slice.\nProceeding with 2D Cellpose segmentation.");
		nucleiSegmentationChoice = "Cellpose segmentation 2D (or on 2D projection)";
	}

	selectWindow(image);
	getDimensions(width, height, channels, slices, frames);
	if(nucleiSegmentationChoice == "Cellpose segmentation 2D (or on 2D projection)") {
		minNucleusSize = PI*Math.sqr((minNucleusSize_setting / pixelWidth / 2));	//Calculate the nucleus area as if it were a circle
		maxNucleusSize = PI*Math.sqr((maxNucleusSize_setting / pixelWidth / 2));	
	}
	else if(nucleiSegmentationChoice == "Cellpose segmentation 3D") {
		minNucleusSize = 4/3*PI*Math.pow((minNucleusSize_setting / pixelWidth / 2), 3);	//Calculate the nucleus area as if it were a sphere
		maxNucleusSize = 4/3*PI*Math.pow((maxNucleusSize_setting / pixelWidth / 2), 3);
	}
	cellpose_input_image = getTitle();
	if(cytoChannel<=0) {
		run("Duplicate...", "title=nuclei duplicate channels=" +nucleiChannel);	//Get nuclei for segmentation
		cellpose_input_image = getTitle();
	}
	else run("Duplicate...", "title=nuclei duplicate channels=" +cytoChannel);	//Get cells (and pretend they are nuclei)

	Cellpose_minCellIntensity = 0;		//To do: maybe make this a parameter for 2D Cellpose as well

	if(slices>1) {
		if(nucleiSegmentationChoice == "Cellpose segmentation 2D (or on 2D projection)") {
			run("Z Project...", "projection=[Max Intensity]");
			rename("forCellpose_Zprojected");
			cellpose_input_image = getTitle();
		}
		else if(nucleiSegmentationChoice == "Cellpose segmentation 3D") {
			if(oldCellposeWrapper == true) exit("Sorry, 3D segmentation with Foci Analyzer is not compatible with the the old BIOP Cellpose wrapper.\nRun Help -> Update... and update the PTBIOP plugins to get the new one.");

			//Get Prefs
			Cellpose_use_GPU = call("ij.Prefs.get", "Cellpose.use.GPU", true);
			Cellpose_3Dseg_type = call("ij.Prefs.get", "Cellpose.3Dseg.type", "2.5D");
			if(CellposeDiameter == 0) CellposeDiameter = call("ij.Prefs.get", "Cellpose.diameter", 30);
			Cellpose_cellprob_threshold = call("ij.Prefs.get", "Cellpose.cellprob.threshold", 0.0);
			Cellpose_anisotropy = call("ij.Prefs.get", "Cellpose.anisotropy", d2s(pixelDepth/pixelWidth,2));
			Cellpose_stitch_threshold = call("ij.Prefs.get", "Cellpose.stitch.threshold", 0.0);
			Cellpose_dP_smooth = call("ij.Prefs.get", "Cellpose.dP.smooth", 0.0);
			Cellpose_min_size = call("ij.Prefs.get", "Cellpose.min.size", 15);
			Cellpose_additional_parameters = call("ij.Prefs.get", "Cellpose.additional.parameters", "");
			Cellpose_normalize3D = call("ij.Prefs.get", "Cellpose.normalize.3D", false);
			Cellpose_minCellIntensity = call("ij.Prefs.get", "Cellpose.min.cell.intensity", 0);

			if(hideCellpose3DDialog == false) {
				Dialog.createNonBlocking("Cellpose 3D advanced settings");
				Dialog.addMessage("Cellpose 3D segmentation parameters", 14, "#000080");
				Dialog.addMessage("If unsure, leave as is or click the Help button below for more info.", 12, "#000080");
				Dialog.addChoice("Segmentation type", newArray("2.5D", "3D"), Cellpose_3Dseg_type);
	//			Dialog.setInsets(-25, 250, 0);
				Dialog.setInsets(-25, 350, 0);
				Dialog.addMessage("'2.5D' means: 2D segmentation on every slice + stitching in z");
				Dialog.addNumber("--diameter", CellposeDiameter, 1, 6, "pixels (N.B. '0 = automatic' does not work for '3D')");
				Dialog.addNumber("--flow_threshold (2.5D only)", probabilityThreshold, 1, 6, "[0-3], default: 0.4");
				Dialog.addNumber("--cellprob_threshold", Cellpose_cellprob_threshold, 1, 6, "[(-6)-6], default: 0.0");
				Dialog.addNumber("--anisotropy", d2s(pixelDepth/pixelWidth,2), 2, 6, "obtained from image dimensions");	//Seems to only work well for 2.5D (?)
				Dialog.addNumber("--stitch_threshold", Cellpose_stitch_threshold, 2, 6, "[0-1], 2D + stitching across planes; only for '2.5D'");
				Dialog.addNumber("--flow3D_smooth (3D only)", Cellpose_dP_smooth, 1, 6, "stddev of Gaussian filter for smoothing flows");
				Dialog.addNumber("--min_size", Cellpose_min_size, 0, 6, "pixels, default: 15");
				Dialog.addString("Additional Cellpose parameters", Cellpose_additional_parameters, 60);
	//			Dialog.setInsets(0, 258, 0);	//probably depends on display scaling and/or font size...
				Dialog.setInsets(0, 290, 0);	//probably depends on display scaling and/or font size...
				Dialog.addCheckbox("--use_gpu", Cellpose_use_GPU);
				Dialog.setInsets(0, 290, 0);
				Dialog.addCheckbox("Normalize the image in 3D before segmentation? (Otherwise possibly leading to falsely detected cells in dim slices.)", Cellpose_normalize3D);
				Dialog.setInsets(10, 0, 0);
				Dialog.addNumber("Remove cells with mean intensity lower than", Cellpose_minCellIntensity, 1, 6, "gray values");
				Dialog.addMessage("Additional Cellpose parameters are e.g. '--no_norm, --restore_type, nuclei, --niter, 200'. Click the Help button for a list of all parameters.");
				Dialog.setInsets(20, 20, 0);
				Dialog.addCheckbox("Hide this dialog for subsequent images", hideCellpose3DDialog);
				Dialog.addHelp("https://cellpose.readthedocs.io/en/latest/cli.html");
				Dialog.show();

				Cellpose_use_GPU = Dialog.getCheckbox();
				Cellpose_3Dseg_type = Dialog.getChoice();
				CellposeDiameter = Dialog.getNumber();
				probabilityThreshold = Dialog.getNumber();
				Cellpose_cellprob_threshold = Dialog.getNumber();
				Cellpose_anisotropy = Dialog.getNumber();
				Cellpose_stitch_threshold = Dialog.getNumber();
				Cellpose_dP_smooth = Dialog.getNumber();
				Cellpose_min_size = Dialog.getNumber();
				Cellpose_additional_parameters = Dialog.getString();
				Cellpose_normalize3D = Dialog.getCheckbox();
				Cellpose_minCellIntensity = Dialog.getNumber();
				hideCellpose3DDialog = Dialog.getCheckbox();
			}

			//Set Prefs
			call("ij.Prefs.set", "Cellpose.use.GPU", Cellpose_use_GPU);
			call("ij.Prefs.set", "Cellpose.3Dseg.type", Cellpose_3Dseg_type);
			call("ij.Prefs.set", "Cellpose.diameter", CellposeDiameter);
			setScriptParameterValue(CellposeDiameter, CellposeDiameter);
			call("ij.Prefs.set", "Cellpose.cellprob.threshold", Cellpose_cellprob_threshold);
			setScriptParameterValue(Cellpose_cellprob_threshold, Cellpose_cellprob_threshold);
			call("ij.Prefs.set", "Cellpose.anisotropy", Cellpose_anisotropy);
			call("ij.Prefs.set", "Cellpose.stitch.threshold", Cellpose_stitch_threshold);
			call("ij.Prefs.set", "Cellpose.dP.smooth", Cellpose_dP_smooth);
			call("ij.Prefs.set", "Cellpose.min.size", Cellpose_min_size);
			call("ij.Prefs.set", "Cellpose.additional.parameters", Cellpose_additional_parameters);
			call("ij.Prefs.set", "Cellpose.normalize.3D", Cellpose_normalize3D);
			call("ij.Prefs.set", "Cellpose.min.cell.intensity", Cellpose_minCellIntensity);
			
			if(Cellpose_use_GPU == true) Cellpose_use_GPU = "--use_gpu";
			if(Cellpose_3Dseg_type == "3D") Cellpose_3Dseg_type = "--do_3D";
			else Cellpose_3Dseg_type = "";
			
			if(Cellpose_normalize3D == true) {
				Ext.CLIJ2_push(cellpose_input_image);
				Ext.CLIJx_normalize(cellpose_input_image, cellpose_input_image_normalized);
				Ext.CLIJ2_pull(cellpose_input_image_normalized);
				Ext.CLIJ2_release(cellpose_input_image);
				Ext.CLIJ2_release(cellpose_input_image_normalized);
				cellpose_input_image = cellpose_input_image_normalized;
				Cellpose_additional_parameters += ", --no_norm";
			}
		}
	}

	setBatchMode("exit and display"); //Cellpose doesn't work in batch mode
	selectWindow(cellpose_input_image);
	showStatus("Performing Cellpose segmentation... (check the Console for live info)");

	//Cellpose segmentation
	//NOTE: for 3D the parameter --dP_smooth has been replaced by --flow3D_smooth (from Cellpose 3.1.something)
	if(nucleiSegmentationChoice == "Cellpose segmentation 2D (or on 2D projection)") {
		if(cytoChannel  >  0) run("Cellpose ...", "env_path="+envPath+" env_type="+envType+" model=["+CellposeModel+"] model_path=["+CellposeModelPath+"] diameter="+CellposeDiameter+" ch1="+cytoChannel+" ch2="+nucleiChannel+" additional_flags=[--use_gpu, --flow_threshold, "+probabilityThreshold+", --cellprob_threshold, 0.0]");
		else run("Cellpose ...", "env_path="+envPath+" env_type="+envType+" model=["+CellposeModel+"] model_path=["+CellposeModelPath+"] diameter="+CellposeDiameter+" ch1="+cytoChannel+" ch2=0 additional_flags=[--use_gpu, --flow_threshold, "+probabilityThreshold+", --cellprob_threshold, 0.0]");
	}
	else if(nucleiSegmentationChoice == "Cellpose segmentation 3D") {
		if(cytoChannel  >  0) run("Cellpose ...", "env_path="+envPath+" env_type="+envType+" model=["+CellposeModel+"] model_path=["+CellposeModelPath+"] diameter="+CellposeDiameter+" ch1="+cytoChannel+" ch2="+nucleiChannel+" additional_flags=["+Cellpose_use_GPU+", "+Cellpose_3Dseg_type+", --flow_threshold, "+probabilityThreshold+", --cellprob_threshold, "+Cellpose_cellprob_threshold+", --stitch_threshold, "+Cellpose_stitch_threshold+", --anisotropy, "+Cellpose_anisotropy+", --flow3D_smooth, "+Cellpose_dP_smooth+", --min_size, "+Cellpose_min_size+", "+Cellpose_additional_parameters+"]");
		else run("Cellpose ...", "env_path="+envPath+" env_type="+envType+" model=["+CellposeModel+"] model_path=["+CellposeModelPath+"] diameter="+CellposeDiameter+" ch1="+nucleiChannel+" ch2=0 additional_flags=["+Cellpose_use_GPU+", "+Cellpose_3Dseg_type+", --flow_threshold, "+probabilityThreshold+", --cellprob_threshold, "+Cellpose_cellprob_threshold+", --stitch_threshold, "+Cellpose_stitch_threshold+", --anisotropy, "+Cellpose_anisotropy+", --flow3D_smooth, "+Cellpose_dP_smooth+", --min_size, "+Cellpose_min_size+", "+Cellpose_additional_parameters+"]");
	}
	else if(oldCellposeWrapper == true) run("Cellpose Advanced", "diameter="+CellposeDiameter+" cellproba_threshold=0 flow_threshold="+probabilityThreshold+" anisotropy=1.0 diam_threshold=12.0 model="+CellposeModel+" nuclei_channel=0 cyto_channel=1 dimensionmode=2D stitch_threshold=-1.0 omni=false cluster=false additional_flags=");

	labelmap_nuclei = getTitle();
	setBatchMode("hide");

//	selectWindow("nuclei");
//	setBatchMode("hide");

	//Close unused images
	if(slices>1 && nucleiSegmentationChoice == "Cellpose segmentation 2D (or on 2D projection)") close("forCellpose_Zprojected");

	Ext.CLIJ2_push(labelmap_nuclei);
	//Exclude nuclei on edges (only XY!) and count nr of nuclei
	if(excludeOnEdges && nucleiSegmentationChoice == "Cellpose segmentation 3D") {
		selectWindow(labelmap_nuclei);
		forEdgeFiltering = "forEdgeFiltering";
		//Make a copy, clear top and bottom slices and get the labels touching the edges
		run("Duplicate...", "title="+forEdgeFiltering+" duplicate");
		Stack.setSlice(1);
		run("Select All");
		run("Clear", "slice");
		Stack.setSlice(slices);
		run("Select All");
		run("Clear", "slice");
		Ext.CLIJ2_push(forEdgeFiltering);
		Ext.CLIJx_flagLabelsOnEdges(forEdgeFiltering, flaglist_labels_on_edges);
		Ext.CLIJ2_release(forEdgeFiltering);
		Ext.CLIJ2_excludeLabels(flaglist_labels_on_edges, labelmap_nuclei, labelmap_nuclei_edges_excluded);
		Ext.CLIJ2_release(flaglist_labels_on_edges);
		Ext.CLIJ2_release(labelmap_nuclei);
	}
	else if(excludeOnEdges && nucleiSegmentationChoice == "Cellpose segmentation 2D (or on 2D projection)") {
		Ext.CLIJ2_push(labelmap_nuclei);
		Ext.CLIJ2_excludeLabelsOnEdges(labelmap_nuclei, labelmap_nuclei_edges_excluded);
		Ext.CLIJ2_release(labelmap_nuclei);
	}
	else if(excludeOnEdges == false) labelmap_nuclei_edges_excluded = labelmap_nuclei;
	close(labelmap_nuclei);

	Ext.CLIJ2_getMaximumOfAllPixels(labelmap_nuclei_edges_excluded, nrNucleiBeforeSizeFiltering);	//Count the number of nuclei
	Ext.CLIJ2_excludeLabelsOutsideSizeRange(labelmap_nuclei_edges_excluded, labelmap_nuclei_filtered, minNucleusSize, maxNucleusSize);
	Ext.CLIJ2_release(labelmap_nuclei_edges_excluded);
	labelmap_nuclei_final = "Labelmap_nuclei";
	Ext.CLIJ2_closeIndexGapsInLabelMap(labelmap_nuclei_filtered, labelmap_nuclei_final);
	Ext.CLIJ2_release(labelmap_nuclei_filtered);
	Ext.CLIJ2_getMaximumOfAllPixels(labelmap_nuclei_final, nrNuclei);	//Count the number of nuclei
	Ext.CLIJ2_pull(labelmap_nuclei_final);
	print(nrNucleiBeforeSizeFiltering - nrNuclei + " cells were removed due to size restrictions (" + minNucleusSize_setting + " - " + maxNucleusSize_setting + " " + unit + " | "+d2s(minNucleusSize,0)+" - "+d2s(maxNucleusSize,0)+" pixels/voxels).");
	if(nrNuclei > 0 && Cellpose_minCellIntensity > 0) {
		if(debugMode) { showImagefromGPU(labelmap_nuclei_final); run("glasbey on dark"); }
		run("Intensity Measurements 2D/3D", "input=nuclei labels=Labelmap_nuclei mean volume");	//Measure intensities using MorphoLibJ
		labels_meanIntensity = Table.getColumn("Mean", "nuclei-intensity-measurements");
		close("nuclei-intensity-measurements");
		//labels_meanIntensity = Table.getColumn("MEAN_INTENSITY", "Results");
		
		labels_meanIntensity = prependToArray(0, labels_meanIntensity);
		Ext.CLIJ2_pushArray(labels_meanIntensity_vector, labels_meanIntensity, labels_meanIntensity.length, 1, 1);
		Ext.CLIJ2_generateParametricImage(labelmap_nuclei_final, labels_meanIntensity_vector, meanIntensity_map);
		Ext.CLIJ2_withinIntensityRange(meanIntensity_map, meanIntensity_map_filtered, Cellpose_minCellIntensity, 1E10);
		Ext.CLIJ2_release(meanIntensity_map);
		//labelmap_nuclei_filtered = "labelmap_nuclei_filtered";
		//labelmap_outlines_filtered = "labelmap_nuclei_filtered_gapsclosed";
		Ext.CLIJ2_multiplyImages(labelmap_nuclei_final, meanIntensity_map_filtered, labelmap_nuclei_filtered);
		Ext.CLIJ2_closeIndexGapsInLabelMap(labelmap_nuclei_filtered, labelmap_nuclei_filtered_gapsclosed);
		Ext.CLIJ2_getMaximumOfAllPixels(labelmap_nuclei_filtered_gapsclosed, nrNuclei);
		Ext.CLIJ2_release(labelmap_nuclei_filtered);
		//Ext.CLIJ2_multiplyImages(nuclei_outlines_org, meanIntensity_map_filtered, nuclei_outlines_filtered);
		//Ext.CLIJ2_pull(nuclei_outlines_filtered);
		//nuclei_outlines = nuclei_outlines_filtered;
		Ext.CLIJ2_release(meanIntensity_map_filtered);
		labelmap_nuclei_final = labelmap_nuclei_filtered_gapsclosed;
		if(debugMode) { showImagefromGPU(labelmap_nuclei_filtered_gapsclosed); run("glasbey on dark"); }
	}
	//Create nuclei outlines from nuclei labelmap
	if (nucleiSegmentationChoice == "Cellpose segmentation 2D (or on 2D projection)") {
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
	}
	else if(nucleiSegmentationChoice == "Cellpose segmentation 3D") {	//Create a 3D image with label edges, without the ugly 'filled planes' that Ext.CLIJ2_detectLabelEdges() outputs.
		Ext.CLIJ2_maximumSliceBySliceSphere(labelmap_nuclei_final, labelmap_maximum, 1, 1);
		Ext.CLIJ2_minimumSliceBySliceSphere(labelmap_nuclei_final, labelmap_minimum, 1, 1);
		Ext.CLIJ2_subtractImages(labelmap_maximum, labelmap_minimum, label_edges);
		Ext.CLIJ2_release(labelmap_maximum);
		Ext.CLIJ2_release(labelmap_minimum);
		Ext.CLIJ2_threshold(label_edges, labelmap_outlines, 1);
		Ext.CLIJ2_pull(labelmap_outlines);
		Ext.CLIJ2_release(label_edges);
		Ext.CLIJ2_release(labelmap_outlines);
		rename("nuclei_outlines");
		run(outlineColor);
	}
	nuclei_outlines = "nuclei_outlines";

	return newArray(labelmap_nuclei_final, nuclei_outlines, nrNuclei);
}


//Manually remove labels by creating point selections. Masks are saved as PNG files (compressed) for later use.
function manually_remove_labels(labelmap, label_edges, nrNuclei, windowName, imageName) {
	selectWindow(windowName);
	Stack.setDisplayMode("composite");
	Stack.setChannel(nucleiChannel);
    setLut(b_reds, b_greens, b_blues);
	Stack.setChannel(fociChannelA);
	run("Green");
   	if(detectFociChannelB) {
   		Stack.setChannel(fociChannelB);
   		run("Red");
   	}
   	run("Add Image...", "image="+label_edges+" x=0 y=0 opacity=100 zero");
	setBatchMode("show");

	if(manualRemoveNuclei == "Load previously saved removals (from output folder)" && File.exists(outputFolder + File.separator + imageName + "__removalMask.png")) {
		print("Previously saved removals found.");
		open(outputFolder + File.separator + imageName + "__removalMask.png");
		if(is("Inverting LUT")) run("Invert LUT");
		run("Points from Mask");
	}
	else if(manualRemoveNuclei == "Load previously saved removals (from output folder)" && !File.exists(outputFolder + File.separator + imageName + "__removalMask.png")) {
		print("[WARNING] No previous removals found. No nuclei will be removed from the analysis.");
	}
	else if(manualRemoveNuclei == "Manually remove nuclei") {
		setTool("multipoint");
		run("Select None");
		run("Point Tool...", "type=Dot color=White size=[Extra Large] show counter=0");
		Overlay.selectable(false)
		waitForUser("Select nuclei to be removed by marking them with a dot by clicking in the image.\n \n* Already placed dots can be shifted by dragging.\n* Alt-click removes a dot.\n* Dots placed outside segmentations have no effect.\n* Press [Shift-a] to start over.\n \nPress OK when finished.");
		if(selectionType == 10) {	//Only create mask if there is a selection
			run("Create Mask");
			saveAs("png", outputFolder + File.separator + imageName+"__removalMask");
			run("Points from Mask");
		}
		else if(File.exists(outputFolder + File.separator + imageName + "__removalMask.png")) {
			File.delete(outputFolder + File.separator + imageName + "__removalMask.png");	//Remove Mask
			print("\\Update:");
		}
	}
	selectWindow(windowName);
	run("Select None");

	close(imageName + "__removalMask");

	Ext.CLIJ2_getMaximumOfAllPixels(labelmap, nrLabels);
	Ext.CLIJ2_pull(labelmap);
	selectWindow(labelmap);
	setTool("rectangle");
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

	nrNuclei_old = nrNuclei;
	Ext.CLIJ2_getMaximumOfAllPixels(labelmap_edited, nrNuclei);
	print(nrNuclei_old - nrNuclei+" were removed manually.");

	label_edges_edited = "label_edges_edited";
	Ext.CLIJ2_detectLabelEdges(labelmap_edited, label_edges_edited);
	Ext.CLIJ2_pull(label_edges_edited);
	run("Multiply...", "value=255");
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

	return newArray(labelmap_edited, label_edges_edited, nrNuclei);
}


function ROI_Manager_to_labelmap(image) {
	setBatchMode(true);
	getDimensions(width, height, channels, slices, frames);
	newImage("labelmap", "16-bit", width, height, 1);
	for (i = 0; i < roiManager("count"); i++) {
		roiManager("select", i);
		changeValues(0, 65535, i+1);
	}
	run("Select None");
	return "labelmap";
}


function detect_foci(image, channel, fociSize, anisotropyFactor, firstTimeProcessing, labelmap_nuclei, labelmap_nuclei_3D, thresholdFactor, threshold, AorB) {
	//Determine threshold on the maximum projection (more robust), only where nuclei are present
	// Small FLAW: the area outside the nuclei is set to zero, but is included in the threshold calculation. (Could be fixed by setting 0 to NaN, but how to do that on the GPU?)
	// On the other hand: there is almost no difference, it seems.
	selectWindow(image);
	run("Select None");
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

	//Determine foci filter sizes
	if(fociSize == "tiny (0.5 pixels)")	{ fociFilterSizeXY = 0.25;	fociSizeXY = 0.5;	}
	if(fociSize == "small (1 pixel)")	{ fociFilterSizeXY = 0.5;	fociSizeXY = 1;		}
	if(fociSize == "average (2 pixels)"){ fociFilterSizeXY = 1.0;	fociSizeXY = 1.5;	}
	if(fociSize == "large (4 pixels)")	{ fociFilterSizeXY = 2;		fociSizeXY = 2;		}
	if(fociSize == "huge (8 pixels)")	{ fociFilterSizeXY = 4;		fociSizeXY = 3;		}
	if(fociSize == "other (define later)") {
		fociSize = getNumber("Enter foci size in pixels", 2);
		fociFilterSizeXY = fociSize/2;	//from diameter ("size") to radius
		if(fociSize <= 1) fociSizeXY = fociSize;
		else if(fociSize <= 2) fociSizeXY = 0.75*fociSize;
		else if(fociSize <= 4) fociSizeXY = 0.66*fociSize;
		else if(fociSize > 4) fociSizeXY = 0.50*fociSize;		
	}

	//pull from GPU, remove outliers, push to GPU
	Ext.CLIJ2_convertFloat(foci, foci_32bit);
	Ext.CLIJ2_pull(foci_32bit);
	foci_outliers_removed = foci+"_outliers_removed";
	rename(foci_outliers_removed);
	selectImage(foci_outliers_removed);
//setBatchMode("show");
	run("Remove Outliers", "block_radius_x="+minOf(25*fociSizeXY, 50)+" block_radius_y="+minOf(25*fociSizeXY, 50)+" standard_deviations=2 stack");
	if(debugMode) setBatchMode("show");
	Ext.CLIJ2_push(foci_outliers_removed);
	if(!debugMode) close(foci_outliers_removed);
//	if(gslices>1) Ext.CLIJ2_mask(foci_outliers_removed, labelmap_nuclei_3D, foci_outliers_removed_masked);
//	else Ext.CLIJ2_mask(foci_outliers_removed, labelmap_nuclei, foci_outliers_removed_masked);
//	showImagefromGPU(foci_outliers_removed_masked);
	//Determine the median of standard deviation signals in the oulier-removed nuclei
	run("Clear Results");
	Ext.CLIJ2_statisticsOfLabelledPixels(foci_outliers_removed, labelmap_nuclei);
	nucleiStddev_ = Table.getColumn("STANDARD_DEVIATION_INTENSITY", "Results");
	medianNucleiStddev = medianOfArray(nucleiStddev_);
	print("Detecting foci in channel "+channel); 
	print("Median value of the Standard Deviations of the nuclear background signal for all nuclei: "+d2s(medianNucleiStddev,1));

	fociFilterSizeZ = 0.5 * anisotropyFactor * fociFilterSizeXY;	//semi-arbitrary, empirically found setting
	fociSizeZ = 0.5 * anisotropyFactor * fociSizeXY;
	//Filter the foci - Difference of Gaussians, but then subtracting a blurred outlier-removed image
	if(gslices>1) Ext.CLIJ2_gaussianBlur3D(foci_32bit, foci_blurred, fociFilterSizeXY/2, fociFilterSizeXY/2, fociFilterSizeZ/2);
	else Ext.CLIJ2_gaussianBlur2D(foci_32bit, foci_blurred, fociFilterSizeXY/2, fociFilterSizeXY/2);
	if(gslices>1) Ext.CLIJ2_gaussianBlur3D(foci_outliers_removed, foci_outliers_removed_blurred, fociFilterSizeXY*4, fociFilterSizeXY*4, fociFilterSizeZ*4);
	else Ext.CLIJ2_gaussianBlur2D(foci_outliers_removed, foci_outliers_removed_blurred, fociFilterSizeXY*4, fociFilterSizeXY*4);

	if(debugMode) showImagefromGPU(foci_blurred);
	if(debugMode) showImagefromGPU(foci_outliers_removed_blurred);
	
	Ext.CLIJ2_release(foci_32bit);
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
	if(processAllOtherImagesFixedThreshold == false) {
		threshold = medianNucleiStddev * thresholdMultiplier;	//Set default threshold at n times the stddev (of the outlier-removed foci, which is a bit lower than the actual stddev, but less dependent on many foci being present)
		threshold = threshold*exp(thresholdFactor);
		print("Automatic threshold set at "+d2s(thresholdMultiplier*exp(thresholdFactor),1)+" times above background stddev (bias factor set at "+thresholdFactor+"): Threshold used: "+d2s(threshold,1));
	}
	else print("Threshold fixed at "+threshold);
	maxDisplaySetting = minOf(pow(2,bits), threshold * 5);

	if(fociDetectionMethod == "AreaMaxima local maximum detection") {
		//Check whether the SCF MPI CBG Plugins is installed
		missingPlugin = "";
		List.setCommands;
		if (List.get("AreaMaxima local maximum detection (2D, 3D)")=="") missingPlugin += "SCF MPI CBG, ";
		if (missingPlugin != "") {
			print("\\Clear");
			missingPlugin = missingPlugin.substring(0, missingPlugin.length-2);
			print("Error: 'AreaMaxima local maximum detection' requires the Fiji Update Site 'SCF MPI CBG' to be activated.\n \nGo to Help -> Update... -> Manage Update Sites and check the relevant box.");
			exit("Error: 'AreaMaxima local maximum detection' requires the Fiji Update Site 'SCF MPI CBG' to be activated.\n \nGo to Help -> Update... -> Manage Update Sites and check the relevant box.\nThis info is also printed to the Log Window.");
		}
		List.clear();

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
		Ext.CLIJ2_mask(labeledSpots, labelmap_nuclei_3D_dilated, labeledSpots_masked);
		if(isOpen("foci_spots_ch"+channel)) close("foci_spots_ch"+channel);
		Ext.CLIJ2_pullBinary(labeledSpots_masked);
		rename("foci_spots_ch"+channel);
		Ext.CLIJ2_release(labeledSpots);
		close(foci_filtered);
		close("Labelmap_detected_foci_ch"+channel);
	}

	else if(fociDetectionMethod == "Marker-controlled watershed (recommended)" || fociDetectionMethod == "Manual thresholding") {
		//	Alternative function - could be faster: 
		//	Ext.CLIJx_detectAndLabelMaximaAboveThreshold(foci, detectedMaxima, 0, 0, 0, threshold, false)
		//	Ext.CLIJ2_pull(detectedMaxima);
		//	rename("detectedMaxima");

		if(fociDetectionMethod == "Manual thresholding") {
			if(!processAllOtherImages && !processAllOtherImagesFixedThreshold) {
				selectImage(image);
				if(firstTimeProcessing == false) setMinAndMax(minDisplayFoci, maxDisplayFoci);
				else {
					getMinAndMax(min, max);
					setMinAndMax(maxOf(0, min), max);
				}
				Ext.CLIJ2_pull(foci_filtered);
				selectImage(foci_filtered);
//				setAutoThreshold("MaxEntropy dark no-reset");	//Unfortunately the threshold window cannot pop up with the 'Don't reset range' button pressed, so this is a bit useless...
				resetMinAndMax;
				getMinAndMax(min, max);
				setMinAndMax(0, max);
				run("16-bit");
				if(firstTimeProcessing == false) setMinAndMax(minDisplayFoci, maxDisplayFoci);
				run("Threshold...");
				image_for_thresholding = foci_filtered + " overlaid with original pixel data (33% opacity)";
				rename(image_for_thresholding);
				if(firstTimeProcessing == true) setThreshold(threshold, pow(2,bits));
				else setThreshold(manualThreshold, pow(2,bits));
				overlay_image_3D(image_for_thresholding, image, 1, 33);
				overlay_image_3D(image_for_thresholding, nuclei_outlines, 1, 50);

				setBatchMode("show");
				waitForUser("Manual thresholding of foci in channel "+channel+":\nThe automatic threshold would be "+d2s(threshold,1)+".\nAdjust the lower threshold value (upper slider) in the Threshold window and press OK to continue");
				if(AorB == "A") {
					getThreshold(manualThresholdA, upper);
					threshold = manualThresholdA;
				}
				else if(AorB == "B") {
					getThreshold(manualThresholdB, upper);
					threshold = manualThresholdB;
				}
				print("\\Update:Overriding threshold with manual value: "+threshold);
				if(optimizationMode == false) processAllOtherImages = true;
			}
			if(processAllOtherImages || processAllOtherImagesFixedThreshold) {
				if(AorB == "A") threshold = manualThresholdA;
				else if(AorB == "B") threshold = manualThresholdB;
				print("\\Update:Using manual threshold set at "+threshold);
			}
		}
		
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
		Ext.CLIJ2_release(maskedSpots);

		Ext.CLIJ2_mask(labeledSpots, labelmap_nuclei_3D_dilated, labeledSpots_masked);
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
			close(foci_filtered);
			close(labeledSpots_masked);
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
		Ext.CLIJ2_release(foci_filtered);
		Ext.CLIJ2_release(foci_filtered_inverted);
//		Ext.CLIJ2_release(labeledSpots_masked);
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
	if(debugMode) showImagefromGPU(labelmapDetectedFociFiltered);

	//Exclude foci outside of nuclei, or nuclei+range
	Ext.CLIJ2_closeIndexGapsInLabelMap(labelmapDetectedFociFiltered, labelmap_foci_final);
	Ext.CLIJ2_release(labelmapDetectedFociFiltered);
	Ext.CLIJ2_getMaximumOfAllPixels(labelmap_foci_final, nrFoci);

	//Mask the spots image with the remaining foci
	Ext.CLIJ2_mask(labeledSpots_masked, labelmap_foci_final, labeledSpots_masked_filtered);
	Ext.CLIJ2_release(labeledSpots_masked);
	if(isOpen("foci_spots_ch"+channel)) close("foci_spots_ch"+channel);
	if(useLargeSpots == true) Ext.CLIJ2_dilateBox(labeledSpots_masked_filtered, labeledSpots_masked_filtered_dilated);
	else labeledSpots_masked_filtered_dilated = labeledSpots_masked_filtered;
	Ext.CLIJ2_pullBinary(labeledSpots_masked_filtered_dilated);
	rename("foci_spots_ch"+channel);

	Ext.CLIJ2_release(labeledSpots_masked_filtered);

	//Remove some lines in the log window (unnecessary info by MorphoLibJ)
	if(fociDetectionMethod == "Marker-controlled watershed (recommended)") {
		logWindowContents = getInfo("log");
		logWindowLines = split(logWindowContents, "\n");
		logWindowLines = Array.trim(logWindowLines, logWindowLines.length - 7);
		logWindowLines = arrayToString(logWindowLines);
		print("\\Clear");
		print(logWindowLines);
	}
	else print("");
	print("\\Update:Channel "+channel+": "+nrFoci+" foci detected in "+nrNuclei+" nuclei ("+d2s(nrFoci/nrNuclei,1)+" foci per nucleus)" + "\n");

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

	return newArray("Labelmap_detected_foci_filtered_ch"+channel, "Mask_foci_ch"+channel, "foci_spots_ch"+channel, threshold);
}


function mergeOriginalAndDetection(original, nrNuclei, nuclei_outlines, foci_mask, foci_spots, channel, useROI, x_ROI, y_ROI) {
	// Merge the original foci channel with the detected foci

	if(useROI == false) {
		//prepare nuclei
		selectWindow("nuclei");
		run(bits+"-bit");	//convert to allow merging
		if(cytoChannel>0 && cellpose==true) {
			rename("cells");
			selectWindow(original);
			run("Duplicate...", "title=nuclei duplicate channels="+nucleiChannel);
			run(bits+"-bit");	//convert to allow merging
		}
		else if(cytoChannel>0) {
			selectWindow(original);
			run("Duplicate...", "title=cells duplicate channels="+cytoChannel);
			run(bits+"-bit");	//convert to allow merging
		}

		//prepare foci raw
		foci_RAW = "foci_ch"+channel;
		if(bits == 16) Ext.CLIJ2_convertUInt16(foci_RAW, foci_RAW_8or16bit);
		else if(bits == 8) Ext.CLIJ2_convertUInt8(foci_RAW, foci_RAW_8or16bit);
		else foci_RAW_8or16bit = foci_RAW;
		Ext.CLIJ2_pull(foci_RAW_8or16bit);
		if(bits == 8 || bits == 16) Ext.CLIJ2_release(foci_RAW_8or16bit);
	}
	if(useROI == true) {
		//prepare nuclei
		close("nuclei");
		selectWindow(original);
		run("Select None");
		run("Duplicate...", "title=nuclei duplicate channels="+nucleiChannel);
		run(bits+"-bit");	//convert to allow merging
		if(cytoChannel>0) {
			selectWindow(original);
			run("Duplicate...", "title=cells duplicate channels="+cytoChannel);
			run(bits+"-bit");
		}

		//prepare foci raw
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
		if(isOpen("Cells")) selectWindow("cells");
		setBatchMode("show");
	}


	if(cytoChannel<=0) run("Merge Channels...", "c1=nuclei c2=" + foci_RAW_8or16bit + " c3=" + foci_mask_large + " c4=" + foci_spots_large + " create keep");
	else run("Merge Channels...", "c1=nuclei c2=" + foci_RAW_8or16bit + " c3=" + foci_mask_large + " c4=" + foci_spots_large + " c5=cells create keep");
	rename("Foci_overlay_ch"+channel);
	getDimensions(width, height, channels, slices, frames);
	Stack.setSlice(slices/2);
	Stack.setChannel(1);
	setLabel("Foci_overlay_ch"+channel, "NUCLEI");
    setLut(b_reds, b_greens, b_blues);
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
	run("Magenta");
	if(cytoChannel>0) {
		Stack.setChannel(5);
		setLabel("Foci_overlay_ch"+channel, "CELLS");
		run("Enhance Contrast...", "saturated=0.35"); 
		setLut(o_reds, o_greens, o_blues);
	}

	//Add nuclei outlines as overlay to the original image
	selectWindow(nuclei_outlines);
	if(overlayBrightness == "dim") setMinAndMax(0, 512);	//Silly way to set nuclei outline overlay to dim, but it works
	if(overlayBrightness == "bright") resetMinAndMax;		//Silly way to set nuclei outline overlay to dim, but it works
	if(debugMode) setBatchMode("show");
	selectWindow("Foci_overlay_ch"+channel);
	for (i = 1; i <= gslices; i++) {
		if(gslices>1 && nucleiSegmentationChoice == "Cellpose segmentation 3D") {
			selectImage(nuclei_outlines);
			Stack.setSlice(i);
			if(maxFociDistanceOutsideNuclei_setting != 0) {
				selectImage(nuclei_dilated_outlines);
				Stack.setSlice(i);
			}
		}
		if(gslices>1) {
			selectImage("Foci_overlay_ch"+channel);
			Stack.setSlice(i);
		}
		run("Add Image...", "image="+nuclei_outlines+" x=0 y=0 opacity="+LABELOPACITY+" zero");	//Add labelmap to image as overlay
		Overlay.setPosition(0, i, 0);
		if(maxFociDistanceOutsideNuclei_setting != 0) run("Add Image...", "image="+nuclei_dilated_outlines+" x=0 y=0 opacity="+maxOf(LABELOPACITY/4, 20)+" zero");
		Overlay.setPosition(0, i, 0);
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
		roiManager("Remove Slice Info");
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


function measureFoci(original, channel, nrNuclei, labelmap_nuclei_3D, labelmap_foci, boundingBox_X, boundingBox_Y, boundingBox_Z, boundingBox_width, boundingBox_height, boundingBox_depth) {
	Ext.CLIJ2_push(labelmap_foci);
//	if(debugMode) showImagefromGPU(labelmap_foci);
	selectWindow(original);
	run("Duplicate...", "title=foci_RAW_ch"+channel + " duplicate channels=" + channel);

//	nucleus_id_ = newArray(nrNuclei);
//	nucleus_area_ =  newArray(nrNuclei);
	foci_count_ = newArray(nrNuclei);
	foci_mean_int_ = newArray(nrNuclei);
	foci_median_int_ = newArray(nrNuclei);
	foci_size_ = newArray(nrNuclei);
	foci_median_size_ = newArray(nrNuclei);
	foci_sum_size_ = newArray(nrNuclei);
	foci_sum_int_ = newArray(nrNuclei);
	nucleus_mean_int_ = newArray(nrNuclei);
	nucleus_sum_intensity_ = newArray(nrNuclei);
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
		labelmap_foci_nucleus_relabeled = "labelmap_foci_nucleus_relabeled";
		foci_raw_cropped = "foci_nucleus";
		// Crop labelmaps and raw images before measuring - much faster for large images
		Ext.CLIJ2_crop3D(foci_raw, foci_raw_cropped, boundingBox_X[i], boundingBox_Y[i], boundingBox_Z[i], boundingBox_width, boundingBox_height, boundingBox_depth);
		if(debugMode) {
//			showImagefromGPU(foci_raw_cropped);
//			getDimensions(width, height, channels, slices, frames);
//			print("dimensions of nucleus "+i+1+" image: "+width+ ", "+height);
		}
		//Crop nucleus (XY only) before measurements, for speed reasons
//		Ext.CLIJ2_crop3D(labelmap_nuclei_3D_dilated, labelmap_nuclei_cropped, boundingBox_X[i], boundingBox_Y[i], 0, boundingBox_width, boundingBox_height, gslices);
//		Ext.CLIJ2_crop3D(labelmap_foci, labelmap_foci_cropped, boundingBox_X[i], boundingBox_Y[i], 0, boundingBox_width, boundingBox_height, gslices);

		Ext.CLIJ2_crop3D(labelmap_nuclei_3D_dilated, labelmap_nuclei_cropped, boundingBox_X[i], boundingBox_Y[i], boundingBox_Z[i], boundingBox_width, boundingBox_height, boundingBox_depth);
		Ext.CLIJ2_crop3D(labelmap_foci, labelmap_foci_cropped, boundingBox_X[i], boundingBox_Y[i], boundingBox_Z[i], boundingBox_width, boundingBox_height, boundingBox_depth);

		Ext.CLIJ2_maskLabel(labelmap_foci_cropped, labelmap_nuclei_cropped, labelmap_foci_nucleus, i+1);

		//This is roughly equally fast as maskLabel:
		//Ext.CLIJ2_labelToMask(labelmap_nuclei_3D, mask_nucleus, i+1);
		//Ext.CLIJ2_mask(labelmap_foci, mask_nucleus, labelmap_foci_nucleus);
		//Also try:
		//Ext.CLIJ2_maskStackWithPlane(Image_source, Image_mask, Image_destination);

		//relabel foci in the current nucleus label and measure using MprpholibJ [faster than Ext.CLIJ2_statisticsOfLabelledPixels(foci_raw_cropped, labelmap_foci_nucleus_relabeled);]
		Ext.CLIJ2_closeIndexGapsInLabelMap(labelmap_foci_nucleus, labelmap_foci_nucleus_relabeled);
		Ext.CLIJ2_pull(foci_raw_cropped);
		Ext.CLIJ2_pull(labelmap_foci_nucleus_relabeled);
		run("Intensity Measurements 2D/3D", "input="+foci_raw_cropped+" labels="+labelmap_foci_nucleus_relabeled+" mean stddev max min median skewness volume");
		fociTable = "foci_nucleus-intensity-measurements";
		fociheadings = split(Table.headings(fociTable), "\t");	//Get column headers. Ok, it is the same for every nucleus, but it's fast anyway.

		//Add intensity measurements, centroid coordinates and circularity/sphericity for all foci
		if(gslices>1) {
			run("Analyze Regions 3D", "sphericity centroid surface_area_method=[Crofton (13 dirs.)] euler_connectivity=6");
			labelTable = "labelmap_foci_nucleus_relabeled-morpho";
		} else {
			run("Analyze Regions", "circularity centroid");
			labelTable = "labelmap_foci_nucleus_relabeled-Morphometry";
		}
		labelheadings = split(Table.headings(labelTable), "\t");	//Get column headers. Ok, it is the same for every nucleus, but it's fast anyway.
		if(Table.size(labelTable) > 0) {
			centroidX_ = Table.getColumn("Centroid.X", labelTable);
			centroidY_ = Table.getColumn("Centroid.Y", labelTable);
			if(gslices>1) centroidZ_ = Table.getColumn("Centroid.Z", labelTable);
			centroidX_ = addScalarToArray(centroidX_, boundingBox_X[i]);
			centroidY_ = addScalarToArray(centroidY_, boundingBox_Y[i]);
			if(gslices>1) addScalarToArray(centroidZ_, boundingBox_Z[i]);
			centroidX_ = multiplyArraywithScalar(centroidX_, pixelWidth);
			centroidY_ = multiplyArraywithScalar(centroidY_, pixelHeight);
			if(gslices>1) centroidZ_ = multiplyArraywithScalar(centroidZ_, pixelDepth);
			Table.setColumn("Centroid.X", centroidX_, labelTable);
			Table.setColumn("Centroid.Y", centroidY_, labelTable);
			if(gslices>1) Table.setColumn("Centroid.Z", centroidZ_, labelTable);
		}
		close(foci_raw_cropped);
		close(labelmap_foci_nucleus_relabeled);
		//Count the number of foci
		Ext.CLIJ2_getMaximumOfAllPixels(labelmap_foci_nucleus_relabeled, nrFoci);
		foci_count_[i] = nrFoci;

		//Copy measured data to allFociResults table
		totalNrFoci = Table.size(allFociResultsTable);
		for (k = totalNrFoci; k < totalNrFoci + nrFoci; k++) {
			Table.set("Cell ID", k, i+1, allFociResultsTable);
			Table.set("channel", k, channel, allFociResultsTable);
		}
	    for (col=0; col<fociheadings.length; col++) {
			col_values = Table.getColumn(fociheadings[col], fociTable);
			for(k = totalNrFoci; k < totalNrFoci + nrFoci; k++) Table.set(fociheadings[col], k, col_values[k - totalNrFoci], allFociResultsTable);
	    }
	    for (col=0; col<labelheadings.length; col++) {
			col_values = Table.getColumn(labelheadings[col], labelTable);
			for(k = totalNrFoci; k < totalNrFoci + nrFoci; k++) Table.set(labelheadings[col], k, col_values[k - totalNrFoci], allFociResultsTable);
	    }

		// calculate stats per nucleus								
		foci_count_[i] = nrFoci;
		if(nrFoci>0) {
			int_ = Table.getColumn("Mean", fociTable);
			foci_mean_int_[i] = meanOfArray(int_);
			foci_median_int_[i] = medianOfArray(int_);
			foci_sum_int_[i] = sumArray(int_);
			size_ = Table.getColumn("Volume", fociTable);
			foci_size_[i] = meanOfArray(size_) * pixelWidth * pixelHeight * pixelDepth;
			foci_median_size_[i] = medianOfArray(size_) * pixelWidth * pixelHeight * pixelDepth;
			foci_sum_size_[i] = foci_size_[i] * nrFoci;	//average size * nr
		}
		else {
			foci_mean_int_[i] = NaN;
			foci_median_int_[i] = NaN;
			foci_sum_int_[i] = NaN;
			foci_size_[i] = NaN;
			foci_median_size_[i] = NaN;
			foci_sum_size_[i] = NaN;
		}
	}
	Table.update(allFociResultsTable);
	
//	//Measure distance to edge of the cell/nucleus
//	//Not active yet: Possibly do this in a loop, because Table.setColumn() overwrites the results of the previous channel.
//	Ext.CLIJ2_distanceMap(labelmap_nuclei_3D, distancemap_nuclei_3D);
//	Ext.CLIJ2_pull(distancemap_nuclei_3D);
//	run("Intensity Measurements 2D/3D", "input="+distancemap_nuclei_3D+" labels="+labelmap_foci+" mean");
//	foci_distance_to_cell_edge_ = Table.getColumn("Mean");
//	Table.setColumn("Distance to cell edge", foci_distance_to_cell_edge_, allFociResultsTable);

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
		nucleus_mean_int_[i] = Table.get("MEAN_INTENSITY", i+1, "Results");
		nucleus_sum_intensity_[i] = Table.get("SUM_INTENSITY", i+1, "Results");
	}
	Ext.CLIJ2_release(foci_raw);

	//Add relevant results to the resultTable
	Table.setColumn("Background intensity ch"+channel, background_, resultTable);
	Table.setColumn("Mean cell intensity ch"+channel, nucleus_mean_int_, resultTable);
	Table.setColumn("Sum cell intensity ch"+channel, nucleus_sum_intensity_, resultTable);
	Table.setColumn("Foci count ch"+channel, foci_count_, resultTable);
	Table.setColumn("Mean foci intensity ch"+channel, foci_mean_int_, resultTable);
	Table.setColumn("Median foci intensity ch"+channel, foci_median_int_, resultTable);
	if(gslices > 1) Table.setColumn("Mean foci volume ("+unit+"^3) ch"+channel, foci_size_, resultTable);
	else Table.setColumn("Mean foci area ("+unit+"^2) ch"+channel, foci_size_, resultTable);
	if(gslices > 1) Table.setColumn("Median foci volume ("+unit+"^3) ch"+channel, foci_median_size_, resultTable);
	else Table.setColumn("Median foci area ("+unit+"^2) ch"+channel, foci_median_size_, resultTable);
	if(gslices > 1) Table.setColumn("Total foci volume ("+unit+"^3) ch"+channel, foci_sum_size_, resultTable);
	else Table.setColumn("Total foci area ("+unit+"^2) ch"+channel, foci_sum_size_, resultTable);
	Table.setColumn("Total foci intensity ch"+channel, foci_sum_int_, resultTable);
	Table.update();
	
	close(fociTable);
	close(labelTable);
	return nrFoci;
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
			Table.setColumn("Mean cell intensity ch"+c, mean_intensity_, resultTable);
			Table.setColumn("Sum cell intensity ch"+c, sum_intensity_, resultTable);
			Table.update;
			Ext.CLIJ2_release(channel_to_measure);
		}
	}
}


//Overlay a 3D image onto another 3D image 
function overlay_image_3D(image, overlay, overlayChannel, opacity) {
	batchMode = is("Batch Mode");
	if(!batchMode) setBatchMode(true);
	selectImage(image);
	getDimensions(width, height, channels, slices, frames);
	Stack.setChannel(overlayChannel);
	Stack.getPosition(channel, slice, frame);
	Overlay.clear;
	for (i = 1; i <= slices; i++) {
		selectImage(overlay);
		Stack.setSlice(i);
		selectImage(image);
		Stack.setSlice(i);
		run("Add Image...", "image=["+overlay+"] x=0 y=0 opacity="+opacity+" zero");
	}
	Stack.setPosition(channel, slice, frame);
	if(!batchMode) setBatchMode(false);
}


//Draw nuclei and/or foci numbers as overlay
function overlay_numbers_on_image(overlay_image) {
	selectImage(overlay_image);
	if(addNumbersOverlay) {
		setFont("SansSerif", labelFontSize, "bold antialiased");
		color1 = color_to_hex(fontColorCells);
		color2 = color_to_hex(fontColorFoci);
		setColor(color1);
		
		if(overlayChoice == "nucleus/cell ID") {
			for (i = 0; i < nrNuclei; i++) {
				x = getResult("CENTROID_X", i);
				y = getResult("CENTROID_Y", i);
				z_start = getResult("BOUNDING_BOX_Z", i);
				z_end = getResult("BOUNDING_BOX_END_Z", i);
				for(k = z_start; k <= z_end ; k++) {
					Overlay.drawString(i+1, x - labelFontSize/2 + x_ROI, y + labelFontSize/2 + y_ROI);
					Overlay.setPosition(0, k+1, 0);
				}
			}
		}
		else if(overlayChoice == "foci count") {
			for (i = 0; i < nrNuclei; i++) {
				x = getResult("CENTROID_X", i);
				y = getResult("CENTROID_Y", i);
				z_start = getResult("BOUNDING_BOX_Z", i);
				z_end = getResult("BOUNDING_BOX_END_Z", i);
				nrf = Table.get("Foci count ch"+fociChannelA, i, resultTable);
				setColor(color2);
				for(k = z_start; k <= z_end ; k++) {
					Overlay.drawString(nrf, x - labelFontSize/2, y + labelFontSize/2);
					Overlay.setPosition(0, k+1, 1);
				}
				if(detectFociChannelB) {
					nrf = Table.get("Foci count ch"+fociChannelB, i, resultTable);
					for(k = z_start; k <= z_end ; k++) {
						Overlay.drawString(nrf, x - labelFontSize/2, y + labelFontSize/2);
						Overlay.setPosition(0, k+1, 2);
					}
				}
			}
		}
		else if(overlayChoice == "both") {
			for (i = 0; i < nrNuclei; i++) {
				x = getResult("CENTROID_X", i);
				y = getResult("CENTROID_Y", i) - labelFontSize/2 - 2;
				z_start = getResult("BOUNDING_BOX_Z", i);
				z_end = getResult("BOUNDING_BOX_END_Z", i);
				setColor(color1);
				for(k = z_start; k <= z_end ; k++) {
					Overlay.drawString(i+1, x - labelFontSize/2, y + labelFontSize/2);
					Overlay.setPosition(0, k+1, 0);
				}
				//y = getResult("CENTROID_Y", i) + labelFontSize/2 + 2;
				y = y + labelFontSize + 4;
				nrf = Table.get("Foci count ch"+fociChannelA, i, resultTable);
				setColor(color2);
				for(k = z_start; k <= z_end ; k++) {
					Overlay.drawString(nrf, x - labelFontSize/2, y + labelFontSize/2);
					Overlay.setPosition(0, k+1, 1);
				}
				if(detectFociChannelB) {
					nrf = Table.get("Foci count ch"+fociChannelB, i, resultTable);
					for(k = z_start; k <= z_end ; k++) {
						Overlay.drawString(nrf, x - labelFontSize/2, y + labelFontSize/2);
						Overlay.setPosition(0, k+1, 2);
					}
				}
			}
		}
		else if(overlayChoice == "none") {
			continue;
		}			
		updateDisplay();
	}
}

function computeOverlap(mask_fociA, mask_fociB, nrNuclei, labelmap_nuclei_3D, boundingBox_X, boundingBox_Y, boundingBox_width, boundingBox_height) {
	showStatus("Computing foci colocalization...");
	foci_overlap_map = "Foci overlap map";
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

	print("\n"+max+" colocalizing foci with overlap > "+ minOverlapSize +" voxels detected in "+nrNuclei+" nuclei ("+d2s(max/nrNuclei,1)+" per nucleus).");

	//Count overlapping foci for all nuclei
	overlap_count_ = newArray(nrNuclei);
	overlap_volume_ = newArray(nrNuclei);
	
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
		Ext.CLIJ2_countNonZeroPixels(labelmap_overlap_nucleus_relabeled);
		overlap_count_[i] = nrOverlappingFoci;
		overlap_volume_[i] = getResult("CountNonZero", i) * pixelWidth * pixelHeight * pixelDepth;

		//Cannot re-use the memory, because the shape is different for every nucleus!
		Ext.CLIJ2_release(labelmap_nuclei_cropped);
		Ext.CLIJ2_release(labelmap_overlap_cropped);
		Ext.CLIJ2_release(labelmap_overlap_nucleus);
		Ext.CLIJ2_release(labelmap_overlap_nucleus_relabeled);
	}
	Ext.CLIJ2_release(labelmap_overlap_filtered);

	//Add overlap results to the resultTable
	Table.setColumn("Overlapping foci count", overlap_count_, resultTable);
	if(gslices > 1) Table.setColumn("Overlapping foci volume ("+unit+"^3)", overlap_volume_, resultTable);
	else Table.setColumn("Overlapping foci area ("+unit+"^2)", overlap_volume_, resultTable);

	//Channel A
	fociCountChannelA_ = Table.getColumn("Foci count ch"+fociChannelA, resultTable);
	overlap_count_fraction_A_ = divideArrays(overlap_count_, fociCountChannelA_);
	overlap_count_percentage_A_ = multiplyArraywithScalar(overlap_count_fraction_A_, 100);
	Table.setColumn("Overlap count % (ch"+fociChannelA+")", overlap_count_percentage_A_, resultTable);

	if(gslices > 1) fociVolumeChannelA_ = Table.getColumn("Total foci volume ("+unit+"^3) ch"+fociChannelA, resultTable);
	else fociVolumeChannelA_ = Table.getColumn("Total foci area ("+unit+"^2) ch"+fociChannelA, resultTable);
	overlap_volume_fraction_A_ = divideArrays(overlap_volume_, fociVolumeChannelA_);
	overlap_volume_percentage_A_ = multiplyArraywithScalar(overlap_volume_fraction_A_, 100);
	if(gslices > 1) Table.setColumn("Overlap volume % (ch"+fociChannelA+")", overlap_volume_percentage_A_, resultTable);
	else Table.setColumn("Overlap area % (ch"+fociChannelA+")", overlap_volume_percentage_A_, resultTable);

	//Channel B
	fociCountChannelB_ = Table.getColumn("Foci count ch"+fociChannelB, resultTable);
	overlap_count_fraction_B_ = divideArrays(overlap_count_, fociCountChannelB_);
	overlap_count_percentage_B_ = multiplyArraywithScalar(overlap_count_fraction_B_, 100);
	Table.setColumn("Overlap count % (ch"+fociChannelB+")", overlap_count_percentage_B_, resultTable);

	if(gslices > 1) fociVolumeChannelB_ = Table.getColumn("Total foci volume ("+unit+"^3) ch"+fociChannelB, resultTable);
	else fociVolumeChannelB_ = Table.getColumn("Total foci area ("+unit+"^2) ch"+fociChannelB, resultTable);
	overlap_volume_fraction_B_ = divideArrays(overlap_volume_, fociVolumeChannelB_);
	overlap_volume_percentage_B_ = multiplyArraywithScalar(overlap_volume_fraction_B_, 100);
	if(gslices > 1) Table.setColumn("Overlap volume % (ch"+fociChannelB+")", overlap_volume_percentage_B_, resultTable);
	else Table.setColumn("Overlap area % (ch"+fociChannelB+")", overlap_volume_percentage_B_, resultTable);


// TO DO (maybe): Use CLIJx_labelOverlapCountMap to calculate the number of overlapping foci.
// It will replace the above code, except for the overlap size filtering, which may not be possible.
// Also, you can use it to count foci in every cell by passing the nuclei labelmap and a foci labelmap.

	//Create Results table column for the overlap visualization, including the background (entry 0)
//	run("Clear Results");
//	setResult("Overlap percentage", 0, 0);
//	for (i = 0; i < nrNuclei; i++) {
//		setResult("Overlap percentage", i+1, overlap_count_percentage_[i]);
//	}
//	updateResults();

	//Visualize the amount of overlapping foci in a map
	//The following line doesn't work. Why not?
//	Ext.CLIJ2_generateParametricImageFromResultsTableColumn(labelmap_nuclei_3D, overlap_count_map, "Overlap percentage");

//	Ext.CLIJ2_pushResultsTableColumn(overlap_fill, "Overlap percentage (");

	overlap_count_map_A = create_overlap_map(overlap_count_percentage_A_, labelmap_nuclei_3D, "overlap_count_percentage_map_A");
	overlap_volume_map_A = create_overlap_map(overlap_volume_percentage_A_, labelmap_nuclei_3D, "overlap_volume_percentage_map_A");
	overlap_count_map_B = create_overlap_map(overlap_count_percentage_B_, labelmap_nuclei_3D, "overlap_count_percentage_map_B");
	overlap_volume_map_B = create_overlap_map(overlap_volume_percentage_B_, labelmap_nuclei_3D, "overlap_volume_percentage_map_B");
	
	if(debugMode) print("Measuring colocalization took "+ getTime() - startTime +" ms.");

	//Create foci mask overlap image 
	run("Merge Channels...", "c1="+mask_fociA+" c2="+mask_fociB+" c3="+overlap_count_map_A+" c4="+overlap_volume_map_A+" c5="+overlap_count_map_B+" c6="+overlap_volume_map_B+" create");
	rename(foci_overlap_map);
	Stack.setChannel(1);
	setLabel(foci_overlap_map, "FOCI ch"+fociChannelA);
	run("Green");
	Stack.setChannel(2);
	setLabel(foci_overlap_map, "FOCI ch"+fociChannelB);
	run("Red");
	Stack.setChannel(3);
	setLabel(foci_overlap_map, "FOCI COUNT OVERLAP PERCENTAGE WITH ch"+fociChannelA);
	setLut(b_reds, b_greens, b_blues);
	setMinAndMax(0, 100);
	Stack.setChannel(4);
	setLabel(foci_overlap_map, "FOCI VOLUME OVERLAP PERCENTAGE WITH ch"+fociChannelA);
	run("Magenta");
	setMinAndMax(0, 100);
	Stack.setChannel(5);
	setLabel(foci_overlap_map, "FOCI COUNT OVERLAP PERCENTAGE WITH ch"+fociChannelB);
	setLut(b_reds, b_greens, b_blues);
	setMinAndMax(0, 100);
	Stack.setChannel(6);
	setLabel(foci_overlap_map, "FOCI VOLUME OVERLAP PERCENTAGE WITH ch"+fociChannelB);
	run("Magenta");
	setMinAndMax(0, 100);
	Stack.setDisplayMode("color");
	setBatchMode("show");

	//Add nuclei outlines as overlay to the overlap image
	selectImage(foci_overlap_map);
	for (c = 1; c <= 6; c++) {
		Stack.setChannel(c);
		for (i = 1; i <= gslices; i++) {
			if(gslices > 1) Stack.setSlice(i);
			run("Add Image...", "image="+nuclei_outlines+" x=0 y=0 opacity="+LABELOPACITY+" zero");	//Add labelmap to image as overlay
		}
	}
	if(addNumbersOverlay) {
		run("Clear Results");
		Ext.CLIJ2_statisticsOfLabelledPixels(labelmap_nuclei, labelmap_nuclei);
		setFont("SansSerif", labelFontSize, "bold antialiased");
		if(overlayChoice == "nucleus/cell ID") {
			color = color_to_hex(fontColorCells);
			setColor(color);
			for (i = 0; i < nrNuclei; i++) {
				x = getResult("MASS_CENTER_X", i);
				y = getResult("MASS_CENTER_Y", i);
				Overlay.drawString(i+1, x - labelFontSize/2, y + labelFontSize/2);
			}
		}
		else if(overlayChoice == "foci count") {
			color = color_to_hex(fontColorCells);
			setColor(color);
			for (i = 0; i < nrNuclei; i++) {
				x = getResult("MASS_CENTER_X", i);
				y = getResult("MASS_CENTER_Y", i);
				Overlay.drawString(Table.get("Foci count ch"+fociChannelA, i, resultTable), x - labelFontSize/2, y + labelFontSize/2);
			}
		}
//TO DO: FIX THIS ALSO FOR BOTH< AND POSSIBLY USE THE FUNCTION overlay_numbers_on_image()
	}
	if(gslices > 1) Stack.setSlice(gslices/2);
	setBatchMode("show");

	return foci_overlap_map;
}


function create_overlap_map(overlap_percentage_, labelmap_nuclei_3D, overlapMapName) {
	overlap_percentage_ = insertElementIntoArrayAtPosition(0, overlap_percentage_, 0);	//Insert element 0
	Ext.CLIJ2_pushArray(overlap_fill, overlap_percentage_, overlap_percentage_.length, 1, 1);
		if(debugMode == true) showImagefromGPU(overlap_fill);
	overlap_map = "overlap_map";
	Ext.CLIJ2_replaceIntensities(labelmap_nuclei_3D, overlap_fill, overlap_map);
	Ext.CLIJ2_release(overlap_fill);
	if(bits == 16) Ext.CLIJ2_convertUInt16(overlap_map, overlap_map_8or16bit);
	else if(bits == 8) Ext.CLIJ2_convertUInt8(overlap_map, overlap_map_8or16bit);
	Ext.CLIJ2_release(overlap_map);
	if(bits == 16) Ext.CLIJ2_replaceIntensity(overlap_map_8or16bit, overlapMapName, 65535, 0);
	else if(bits == 8) Ext.CLIJ2_replaceIntensity(overlap_map_8or16bit, overlapMapName, 255, 0);
	Ext.CLIJ2_pull(overlapMapName);
	rename(overlapMapName);
	if(debugMode == true) { run("Duplicate...", "title="+overlapMapName+"_debug"); setBatchMode("show"); }
	Ext.CLIJ2_release(overlapMapName);
	Ext.CLIJ2_release(overlap_map_8or16bit);
	
	return overlapMapName;
}	


//Create a nice azure blue LUT for the nuclei
function create_azure_lut(b_reds, b_greens, b_blues) {
	for (i=0; i<256; i++) {
		b_reds[i] = 0;
		b_greens[i] = 0.5*i;
		b_blues[i] = i;
	}
}

//Create an orange LUT for the cells
function create_orange_lut(o_reds, o_greens, o_blues) {
	for (i=0; i<256; i++) {
		o_reds[i] = i;
		o_greens[i] = 0.5*i;
		o_blues[i] = 0;
	}
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


//Prepends the value to the array
function prependToArray(value, array) {
	temparray=newArray(lengthOf(array)+1);
	for (i=0; i<lengthOf(array); i++) {
		temparray[i+1]=array[i];
	}
	temparray[0]=value;
	array=temparray;
	return array;
}


//Adds a scalar to all elements of an array
function addScalarToArray(array, scalar) {
	added_array=newArray(lengthOf(array));
	for (a=0; a<lengthOf(array); a++) {
		added_array[a]=array[a] + scalar;
	}
	return added_array;
}


//Returns the sum of all elements of an arrays, ignoring NaNs
function sumArray(array) {
	sum=0;
	for (a=0; a<lengthOf(array); a++) {
		if(!isNaN(array[a])) sum=sum+array[a];
	}
	return sum;
}


//Returns the mean of the array
function meanOfArray(array) {
	Array.getStatistics(array, min, max, mean, stdDev);
	return mean;
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


//Insert an element with value at a certain position
function insertElementIntoArrayAtPosition(value, array, position) {
	if (position<lengthOf(array)) {
		Array.rotate(array, -position);
		Array.reverse(array);
		array[array.length]=value;
		Array.reverse(array);
		Array.rotate(array, position);
	}
	else array[array.length]=value;
	return array;
}


function showImagefromGPU(string_Image) {
	Ext.CLIJ2_pull(string_Image);
	setBatchMode("show");
}


//Get the persistent value of the script parameter 'param' in class. N.B. This returns 'null' when the parameter is set to the default value!
function getPref(class, param) {
	return eval("js",
		"var ctx = Packages.ij.IJ.runPlugIn(\"org.scijava.Context\", \"\");" +
		"var ps = ctx.service(Packages.org.scijava.prefs.PrefService.class);" +
		"var " + param + " = ps.get(" + class + ".class, \"" + param + "\", \"<null>\");" +
		param + ";"
	);
}


//Gets a persistent Script Parameter value. N.B. This returns 'null' when the parameter is set to the default value!
function getScriptParameterValue(string_parameter) {
	return eval("js",
		"var ctx = Packages.ij.IJ.runPlugIn(\"org.scijava.Context\", \"\");" +
		"var ps = ctx.service(Packages.org.scijava.prefs.PrefService.class);" +
		"var " + string_parameter + " = ps.get(org.scijava.script.ScriptModule.class, \"" + string_parameter + "\", \"<null>\");" +
		string_parameter + ";"
	);
}


//Sets a persistent Script Parameter value. N.B. This returns 'null' when the parameter is set to the default value!
function setScriptParameterValue(string_parameter, string_value) {
	eval("js",
		"var ctx = Packages.ij.IJ.runPlugIn(\"org.scijava.Context\", \"\");" +
		"var ps = ctx.service(Packages.org.scijava.prefs.PrefService.class);" +
		"var " + string_parameter + " = ps.put(org.scijava.script.ScriptModule.class, \"" + string_parameter + "\", \"" + string_value + "\");"
	);
}


function cleanup() {
//	close("nuclei");
	close("foci");
}

