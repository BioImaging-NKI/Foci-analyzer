import qupath.lib.common.ColorTools
import qupath.lib.gui.scripting.QPEx
import qupath.lib.objects.classes.PathClass
import qupath.lib.scripting.QP
import qupath.lib.classifiers.object.ObjectClassifiers.ClassifyByMeasurementBuilder
import qupath.opencv.ml.pixel.PixelClassifierTools

//user parameters
def MEASUREMENT_NAME = "Foci count ch3"  // for positive/negative classification
def MIN_THRESHOLD = 2                    // for positive/negative classification
def POSITIVE_COLOR = ColorTools.packRGB(255,0,0)
def NEGATIVE_COLOR = ColorTools.packRGB(0,255,0)
def CLEAR_MEASUREMENTS = true            // Clear all measurements
def TILE_IDENTIFIER = "exported_tiles"   // Annotations with this name will be ignored
def RUN_ON = RunOn.ALL_ANNOTATIONS       // Other options: RunOn.SELECTED or RunOn.ALL_ANNOTATIONS or RunOn.WHOLE_IMAGE
def SEARCH_HEADER = 'foci'               // checks the .tsv header if it contains this character sequence, and ignores case.
                                         // these are added as measurements to the dectections.

def mydets = QP.getDetectionObjects()

if (CLEAR_MEASUREMENTS){mydets.parallelStream().forEach{it.measurementList.clear()}}

// get files per annotation
FileFilter filter = new FileFilter() {
    public boolean accept(File f)
    {
        return f.getName().endsWith("_Foci_results.tsv")
    }
}
def pathObjects = QP.getSelectedObjects()
if (RUN_ON == RunOn.ALL_ANNOTATIONS || RUN_ON == RunOn.WHOLE_IMAGE) {
    pathObjects = QP.getAnnotationObjects()
            .findAll { it.hasROI() & it.getROI().isArea() }
            .findAll { !(it.getName() == TILE_IDENTIFIER) }
}
def outpath = QP.buildFilePath(QP.PROJECT_BASE_DIR, "Exported_tiles", QP.getCurrentImageNameWithoutExtension())
pathObjects.each { anno ->
    def curr_inpath = new File(QP.buildFilePath(outpath, anno.getID().toString()))
    if (!curr_inpath.exists()){
        print("ERROR: "+curr_inpath.toString()+" not found.")
        return
    }
    def tsvfiles = curr_inpath.listFiles(filter)
    print("Found "+tsvfiles.size()+ " tsv-files.")
    tsvfiles.each {file->
        if (!file.exists()) {
            print("File not Found " + file.toString())
            return
        }
        def data = file.readLines()*.split('\t')
        def header = data[0]
        def measures_to_add = header.findAll {it.containsIgnoreCase(SEARCH_HEADER)}
        def uuid_idx = header.findIndexOf {it=="Cell UUID"}
        def idxs_to_add = measures_to_add.collect {measure-> header.findIndexOf {it==measure}}

        // get pathobjects
        def allobjects = anno.getChildObjects()
        data.eachWithIndex{ String[] tsvline, int tsvlineidx ->
            if (tsvlineidx>0){
                def uuid = UUID.fromString(tsvline[uuid_idx])
                def found_object = allobjects.find({it.getID()==uuid})
                if (found_object==null){
                    print("WARNING: not found "+uuid.toString())
                }else{
                    def measlist = found_object.getMeasurementList()
                    idxs_to_add.each {
                        measlist.put(header[it], Double.parseDouble(tsvline[it]))
                    }
                    measlist.close()
                }
            }
        }
    }
}

QP.fireHierarchyUpdate()

// Classify nuclei on a foci measurement
def classifier = new ClassifyByMeasurementBuilder(MEASUREMENT_NAME)
    .below("Negative")
    .aboveEquals("Positive")
    .threshold(MIN_THRESHOLD)
    .build()
classifier.classifyObjects(QP.getCurrentImageData(),true)

def project = QP.getProject()
def positive = project.getPathClasses().find {it.getName()=="Positive"}
if (positive==null){
    def newclasses = project.getPathClasses().collect()
    newclasses.add(PathClass.fromString("Positive"))
    project.setPathClasses(newclasses)
    project.syncChanges()
    positive = project.getPathClasses().find {it.getName()=="Positive"}
}
positive.setColor(POSITIVE_COLOR)
def negative = project.getPathClasses().find {it.getName()=="Negative"}
if (negative==null){
    def newclasses = project.getPathClasses().collect()
    newclasses.add(PathClass.fromString("Negative"))
    project.setPathClasses(newclasses)
    project.syncChanges()
    negative = project.getPathClasses().find {it.getName()=="Negative"}
}
negative.setColor(NEGATIVE_COLOR)
QPEx.getQuPath().refreshProject()
print(project.getPathClasses())

def detobjs = QP.getDetectionObjects()
int nr_foci = 0
int nr_pos_cells = 0
detobjs.each {
    def measlist = it.getMeasurementList()
    nr_foci += (int) measlist.get(MEASUREMENT_NAME)
    if (it.getPathClass().toString()=="Positive"){nr_pos_cells++}
}
print("nr of nuclei: "+detobjs.size())
print("nr of positive nuclei: "+nr_pos_cells)
print("nr of foci: "+nr_foci)
print("avg nr of foci per nucleus: "+nr_foci/detobjs.size())

// CLASS DEFINITIONS
enum RunOn {
    SELECTED, ALL_ANNOTATIONS, WHOLE_IMAGE
}
