import qupath.lib.images.servers.ImageServer
import qupath.lib.objects.PathObject
import qupath.lib.objects.PathObjectTools
import qupath.lib.objects.classes.PathClass
import qupath.lib.objects.classes.Reclassifier
import qupath.lib.roi.interfaces.ROI
import qupath.lib.scripting.QP
import qupath.opencv.ml.pixel.PixelClassifierTools

import java.awt.image.BufferedImage

def PIXEL_CLASSIFIER = "hypoxia"
def imageData = QP.getCurrentImageData()
def classifier = QP.getProject().getPixelClassifiers().get(PIXEL_CLASSIFIER)
def classifierServer = PixelClassifierTools.createPixelClassificationServer(imageData, classifier)
classifyObjectsByCentroid(classifierServer, QP.getDetectionObjects(), true)

public static void classifyObjectsByCentroid(ImageServer<BufferedImage> classifierServer, Collection<PathObject> pathObjects, boolean preferNucleusROI) {
    Map<Integer, PathClass> labels = classifierServer.getMetadata().getClassificationLabels();
    List<Reclassifier> reclassifiers = pathObjects.parallelStream().map((p) -> {
        try {
            ROI roi = PathObjectTools.getROI(p, preferNucleusROI);
            int x = (int)roi.getCentroidX();
            int y = (int)roi.getCentroidY();
            int ind = PixelClassifierTools.getClassification(classifierServer, x, y, roi.getZ(), roi.getT());
            return new Reclassifier(p, (PathClass)labels.getOrDefault(ind, (Object)null), true);
        } catch (Exception var8) {
            return new Reclassifier(p, (PathClass)null, true);
        }
    }).toList();
    reclassifiers.parallelStream().forEach((r) -> {
        r.apply();
    });
}
print("Done running pixel classifier: "+PIXEL_CLASSIFIER)