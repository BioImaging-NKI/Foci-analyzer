import qupath.ext.stardist.StarDist2D
import qupath.lib.regions.RegionRequest
import qupath.lib.roi.ROIs
import qupath.lib.images.writers.TileExporter
import com.google.gson.Gson
import com.google.gson.stream.JsonReader
import qupath.lib.objects.PathObjects
import qupath.lib.objects.PathObject
import ij.plugin.frame.RoiManager
import qupath.imagej.tools.IJTools
import qupath.lib.gui.QuPathGUI
import qupath.lib.scripting.QP


//user parameters
def RUN_ON = RunOn.ALL_ANNOTATIONS          // Other options: RunOn.SELECTED or RunOn.ALL_ANNOTATIONS or RunOn.WHOLE_IMAGE
def MODELPATH = 'D:\\Software\\StarDist_pretrained_models\\dsb2018_heavy_augment.pb'
def CHANNEL = 'Channel 1'                   // channel name for StarDist segmentation
def MIN_NUCLEUS_SIZE = 20.0                 // um^2
def TILE_SIZE = 1600                        // in pixels
def OVERLAP = 128                           // tile overlap in pixels
def INCLUDE_PARTIAL_TILES = true            // Specify whether incomplete tiles at image boundaries should be included
def ADD_TILES_AS_ANNOTATIONS = true         // Show exported tiles as annotations. Existing tile annotations will be removed upon running the script.
def REMOVE_ALL = true                       // remove all detections at the start
def TILE_IDENTIFIER = "exported_tiles"      // Annotations with this name will be removed at the start

def gui = QuPathGUI.getInstance()
def options = gui.getOverlayOptions()
options.setShowNames(false)
// options.setOpacity(0.7)

if (REMOVE_ALL){
    QP.removeDetections()
}
if (RUN_ON == RunOn.WHOLE_IMAGE){
    QP.createFullImageAnnotation(true)
    RUN_ON = RunOn.SELECTED
}
// remove the possible earlier rectangles created by this script
def removerects = QP.getAnnotationObjects().findAll {it.getName()==TILE_IDENTIFIER}
QP.removeObjectsAndDescendants(removerects)

// Customize how the StarDist detection should be applied.
// This uses a 'builder' to make it easier to add lots of options.
// IMPORTANT! You probably don't need all these -
// read the descriptions & remove the lines you don't want

def stardist = StarDist2D.builder(MODELPATH)
        .normalizePercentiles(1, 99)
        .channels(CHANNEL)                 // Select detection channel (usually useful for fluorescence, not needed for RGB);
                                           // the channel can be selected by name or index/number (where 0 is the first channel)
        .threshold(0.5)                    // Probability (detection) threshold
        .pixelSize(0.5)                    // Resolution for detection
        .tileSize(1024)                    // Specify width & height of the tile used for prediction
//        .cellExpansion(5.0)              // Approximate cells based upon nucleus expansion
        .cellConstrainScale(1.5)           // Constrain cell expansion using nucleus size
        .ignoreCellOverlaps(false)         // Set to true if you don't care if cells expand into one another
//        .measureShape()                  // Add shape measurements
        .measureIntensity()                // Add cell measurements (in all compartments)
//        .includeProbability(true)        // Add probability as a measurement (enables later filtering)
//        .nThreads(4)                     // Limit the number of threads used for (possibly parallel) processing
//        .simplify(1)                     // Control how polygons are 'simplified' to remove unnecessary vertices
//        .doLog()                         // Use this to log a bit more information while running the script
//        .createAnnotations()             // Generate annotation objects using StarDist, rather than detection objects
//        .constrainToParent(false)        // Prevent nuclei/cells expanding beyond any parent annotations (default is true)
//        .classify("Tumor")               // Automatically classify all created objects as 'Tumor'
        .build()

// Define which objects will be used as the 'parents' for detection
def pathObjects = QP.getAnnotationObjects().findAll {it.hasROI() & it.getROI().isArea()}
if (RUN_ON==RunOn.SELECTED) {pathObjects = QP.getSelectedObjects()}

// Run detection for the objects
print("Running Stardist on "+pathObjects.size()+" annotations.")
def imageData = QP.getCurrentImageData()
if (pathObjects.isEmpty()) {
    QP.getLogger().error("No annotations found to run analysis on!")
    return
}
stardist.detectObjects(imageData, pathObjects)
stardist.close() // This can help clean up & regain memory

// Remove nuclei that are too small
def nuclei = QP.getDetectionObjects()
def pixelSize = QP.getCurrentImageData().getServer().getPixelCalibration().getAveragedPixelSize()
def remove = nuclei.findAll({
    !(it.getROI().getScaledArea(pixelSize as double, pixelSize as double) > MIN_NUCLEUS_SIZE)
})
QP.removeObjectsAndDescendants(remove)
QP.fireHierarchyUpdate()

def outpath = QP.buildFilePath(QP.PROJECT_BASE_DIR, "Exported_tiles", QP.getCurrentImageNameWithoutExtension())
pathObjects.each {anno->
    def curr_outpath = new File(QP.buildFilePath(outpath,anno.getID().toString()))
    curr_outpath.mkdirs()
    def region = RegionRequest.createInstance(QP.getCurrentServer() as String, 1.0, anno.getROI())
    new TileExporter(QP.getCurrentImageData())
            .overlap(OVERLAP)  // pixels
            .tileSize(TILE_SIZE)
            .region(region)
            .annotatedTilesOnly(true)
            .includePartialTiles(INCLUDE_PARTIAL_TILES)
            .exportJson(true)
            .writeTiles(curr_outpath as String)
}
def pixelsize = QP.getCurrentServer().getPixelCalibration().getAveragedPixelSize().doubleValue()
def pixelsizepath = new File(QP.buildFilePath(outpath, 'pixelsize.txt'))
pixelsizepath.write(pixelsize.toString())

//Load JSON files containing tile locations and save detections as imagej roi
double downsample = 1.0
pathObjects.each {anno->
    def curr_outpath = QP.buildFilePath(outpath,anno.getID().toString())
    def jsonpath = new File(QP.buildFilePath(curr_outpath, QP.getCurrentImageNameWithoutExtension()+'-tiles.json'))
    JsonReader reader = new JsonReader(new FileReader(jsonpath))
    Gson gson = new Gson()
    Tiles tiles = gson.fromJson(reader, Tiles)
    reader.close()

    //Add tiles from JSON
    Collection<PathObject> rectangles = []
    tiles.tiles.forEach {
        def roi = ROIs.createRectangleROI(it.region.getX() as double, it.region.getY() as double, it.region.getWidth() as double, it.region.getHeight() as double, anno.getROI().getImagePlane())
        def rect = PathObjects.createAnnotationObject(roi)
        rect.setName(TILE_IDENTIFIER)
        rectangles.add(rect)
    }
    ArrayList<Collection<PathObject>> detection_reference = new ArrayList<Collection<PathObject>>()
    for (int i=0 ; i<rectangles.size(); i++){
        detection_reference[i] = new ArrayList<PathObject>()
    }
    print("Processing "+QP.getDetectionObjects().size()+" detections")
    QP.getDetectionObjects().eachWithIndex { PathObject det, int i ->
        def idx = rectangles.findIndexValues {it.getROI().contains(det.getROI().getCentroidX(), det.getROI().getCentroidY())}
        idx.each {
            detection_reference[it].add(det)
        }
    }
    def roiMan = new RoiManager(false)
    rectangles.eachWithIndex { PathObject rectangle, int i ->
        print("Processing tile "+(i+1)+"/"+rectangles.size())
        def roi = rectangle.getROI()
        def name = QP.getCurrentImageNameWithoutExtension() + " [x="+roi.getBoundsX().intValue()+",y="+roi.getBoundsY().intValue()+",w="+roi.getBoundsWidth().intValue()+",h="+roi.getBoundsHeight().intValue()+"]_ROIs.zip"
        double x = -roi.getBoundsX()
        double y = -roi.getBoundsY()
        detection_reference[i].each {det->
            def IJroi = IJTools.convertToIJRoi(det.getROI(), x, y, downsample)
            IJroi.setName(det.getID().toString())
            roiMan.addRoi(IJroi)
        }
        def path = QP.buildFilePath(curr_outpath, name)
        if(detection_reference[i].size() > 0) roiMan.runCommand("Save", path)
        roiMan.reset()
    }
    if(ADD_TILES_AS_ANNOTATIONS) QP.addObjects(rectangles)
    rectangles.each {it.setColor(255,255,255)}

    print("Tiles exported to "+outpath)
}
println('Done!')


// CLASS DEFINITIONS
class Tiles {
    String qupath_version
    String base_directory
    List<Tile> tiles
}

class Tile {
    Region region
    String image
}

class Region {
    String path
    double downsample
    int x
    int y
    int width
    int height
    int z
    int t
}

enum RunOn {
    SELECTED, ALL_ANNOTATIONS, WHOLE_IMAGE
}
