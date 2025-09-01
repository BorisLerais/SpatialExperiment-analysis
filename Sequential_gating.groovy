// The general idea is to apply the combined classifier sequentially
// starting from the cleanest classifier and attribute class we are most confident with
// After that we remove the classified cells from the detections list and proceed to the next classifier
// This way we can keep the classification of B cells while working on T cells and so on

// START-UP
// get the image data and hierarchy
//detections is the list of Detection objects, here our cells
import qupath.lib.gui.scripting.QPEx
def imageData = getCurrentImageData()
def hierarchy = imageData.getHierarchy()
def detections = hierarchy.getDetectionObjects()

// We define three methods : 
// removeClass is a helper method that takes a detectionList and a className
// Then return the detectionList after removing all cells of the given class

// gateDetections is the main method that takes a detectionList, a list of phenotype and a class name
// it first attributes the class to each cell corresponding to one of the phenotypes
// We also bind the cells to global, to be added to the filterList
// then removes the cells from the detectionList

// nexClassifier is a work around to only reclassify non filtered cells
// it records the classification of the already filtered cells
// then apply a classifier
// and finally restores the previous class of filtered cells
def removeClass(detectionList, className) { 
    detectionList.removeAll { detection ->
        detection.getPathClass()?.getName() == className
    }
    return detectionList  // for method chaining if desired
    }


def gateDetections(detectionList, phenoList, className) {
    def cells = detectionList.findAll { obj ->
            obj.getClassification() in phenoList
        }
    def pathClass = getPathClass(className)
    
    // We set a classification for the cells of corresponding phenotype
    cells.forEach { it.setPathClass(pathClass)}
    
    binding.setVariable(className, cells)
    
    // And we remove the cells from the collection
    removeClass(detectionList, className)
    
    return detectionList
    }

def nextClassifier(classifierName, filteredList) {
    //Remember old class of previously classified objects
    def oldClass = [:]
    
    filteredList.each { obj ->
        oldClass[obj] = obj.getPathClass()
        }
    // Run the next classifier
    runObjectClassifier(classifierName)
    
    // Restore previous classification
    filteredList.each { obj ->
        obj.setPathClass(oldClass[obj])
        }
    }

// filteredList will contain the cells that were already classified
def filteredList = []

// B CELLS GATE :
// Gating based on CD45, CD3 and CD20
// We derive B_cells and BnT_cells class from it (T cells are gated later as CD4 or CD8 T cells)
runObjectClassifier("B_cells")

// Define the phenotypes : 
def B_pheno = [
    "CD45pos: CD3neg: CD20pos"
    ]
def BnT_pheno = ["CD45pos: CD3pos: CD20pos"]

// Gate the cells and add them to the filteredList
gateDetections(detections, B_pheno, "B_cells")
filteredList.addAll(B_cells)
gateDetections(detections, BnT_pheno, "BnT_cells")
filteredList.addAll(BnT_cells)

// T CELLS GATE :
// Gating based on CD45, CD3, CD20, CD4 and CD8
// We derive CD4_T_cells and CD8_T_cells from this gate, as well as Double positive (DP) and Double negative (DN)
// Double positive T cells might represent a group similar to BnT with interaction of a CD4 and a CD8 cell
nextClassifier("T_cells", filteredList)

def CD4_T_pheno = [
        "CD45pos: CD3pos: CD20neg: CD4 high: CD8 neg",
        "CD45pos: CD3pos: CD20neg: CD4 low: CD8 neg"
    ]
def CD8_T_pheno = [
        "CD45pos: CD3pos: CD20neg: CD4 neg: CD8 pos"
    ]
def DP_T_pheno = [
        "CD45pos: CD3pos: CD20neg: CD4 high: CD8 pos",
        "CD45pos: CD3pos: CD20neg: CD4 low: CD8 pos"
    ]

gateDetections(detections, CD8_T_pheno, "CD8_T_cells")
filteredList.addAll(CD8_T_cells)
gateDetections(detections, CD4_T_pheno, "CD4_T_cells")
filteredList.addAll(CD4_T_cells)
gateDetections(detections, DP_T_pheno, "DP_T_cells")
filteredList.addAll(DP_T_cells)

// Macrophages
// Gating based on CD45, CD3, CD20, CD4 and CD68
// We only derive Macrophages from this gates
// Macrophages express low CD4 level, we will accept both CD4neg and CD4low, but not CD4 high and never CD3pos
nextClassifier("Macrophages", filteredList)

def macro_pheno = [
    "CD45pos: CD3neg: CD20neg: CD4 low: CD68pos",
    "CD45pos: CD3neg: CD20neg: CD4 neg: CD68pos"
    ]

gateDetections(detections, macro_pheno, "Macrophages")
filteredList.addAll(Macrophages)

// Plasma cells
// Gating based of CD45, CD3, CD20, and CD138
// We derive plasma_cells and CD138low from this gate
nextClassifier("Plasma_cells", filteredList)

def plasma_pheno = [
    "CD45pos: CD3neg: CD20neg: CD138 high"
    ]
// CD138low cells can include endothelial cells, epithelial cells, even fibroblast etc.
def CD138_low_pheno = [
    "CD45neg: CD3neg: CD20neg: CD138 low"
    ]

gateDetections(detections, plasma_pheno, "Plasma_cells")
filteredList.addAll(Plasma_cells)
gateDetections(detections, CD138_low_pheno, "CD138low")
filteredList.addAll(CD138low)

// NK cells
// Gating based on CD45, CD3, CD20,  CD68, CD138 and CD56, only deriving CD56+ cells as NK

nextClassifier("NK_cells", filteredList)
def NK_pheno = [
        "CD45pos: CD3neg: CD20neg: CD68neg: CD138 neg: CD56pos"
    ]

gateDetections(detections, NK_pheno, "NK_cells")
filteredList.addAll(NK_cells)

// Synoviocytes :
// The rest of immune cells without specific marker will be classified in a other_immune class (CD4low and CD138low accepted)
// Cells with no specific marker expression will be classified as synoviocytes. CD4low will be accepted due to low specificity
// CD138low with no other marker has been addressed with plasma cells
nextClassifier("Main_populations", filteredList)

def immune_pheno = [
    "CD45pos: CD20neg: CD3neg: CD8 neg: CD4 neg: CD68neg: CD138 neg: CD56neg",
    "CD45pos: CD20neg: CD3neg: CD8 neg: CD4 low: CD68neg: CD138 neg: CD56neg",
    "CD45pos: CD20neg: CD3neg: CD8 neg: CD4 neg: CD68neg: CD138 low: CD56neg",
    "CD45pos: CD20neg: CD3neg: CD8 neg: CD4 low: CD68neg: CD138 low: CD56neg"
    ]

def syno_pheno = [
    "CD45neg: CD20neg: CD3neg: CD8 neg: CD4 neg: CD68neg: CD138 neg: CD56neg",
    "CD45neg: CD20neg: CD3neg: CD8 neg: CD4 low: CD68neg: CD138 neg: CD56neg"
    ]
gateDetections(detections, immune_pheno, "other_immune")
filteredList.addAll(other_immune)
gateDetections(detections, syno_pheno, "Synoviocytes")
filteredList.addAll(Synoviocytes)

// rest will be classified as unknown
// another script will try to adress these unknown phenotype to their most likely class
def Unk_class = getPathClass("Unknown")
detections.forEach {it.setPathClass(Unk_class)}

// Now we will store the class of all cells
// So we will take all the detections again :
detections = hierarchy.getDetectionObjects()

// Store the class as a csv :
mkdirs(buildFilePath(PROJECT_BASE_DIR, "exports"))
mkdirs(buildFilePath(PROJECT_BASE_DIR, "exports/sequential_gating"))
mkdirs(buildFilePath(PROJECT_BASE_DIR, "exports/sequential_gating/minimal"))
def outputFile = buildFilePath(PROJECT_BASE_DIR, "exports/sequential_gating/minimal", "${QPEx.getProjectEntry().getImageName()}_cells.csv")

// Open file for writing
def writer = new File(outputFile)
writer.withPrintWriter { pw ->
    // Write header
    pw.print("Name,Class")
    pw.println()

    // Loop through detections (cells)
    detections.eachWithIndex { det, i ->
        def name = "Cell_${i+1}"
        def pathClass = det.getPathClass()?.getName() ?: "None"

        pw.print("${name},${pathClass}")
        pw.println()
    }
}

print "Exported ${detections.size()} cells to: ${outputFile}"