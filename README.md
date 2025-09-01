# Anatomy-IntersectionDistribution-Analysis

# MATLAB classes for analyzing intersections and distributions of 3D defect models

Developed in the context of orthopaedic research for pelvic defect evaluation, but also applicable to other anatomical or technical surface models.

- **Intersection.m** – computation of intersections between defect volume models, intersection thresholds  
- **Distribution.m** – statistical evaluation and visualization (frequencies)

---

## Overview
MATLAB classes for analyzing and visualizing intersections. Includes tools for 3D overlap computation and interactive slice/3D distribution views.

### Key Features

#### Intersection
1. **Sequential intersections:** `intersectSeq` – Computes sequential intersections across combinations (not recommended)  
2. **Pairwise Intersections** for any intersection combinations (recommended)
   - `intersectBox`: Pre-filters points using a non-axis-aligned bounding box  
   - `intersectPairs`: Performs pairwise point-in-volume tests using *inpolyhedron*, returns index masks  
   - `intersectCount`: Builds frequency tables across pairings and filters by intersection thresholds  
3. **Comparison:** Compares intersection types (unique vs. shared vertices/points)  

#### Distribution
1. **Distribution:** Distribution of anatomical models with reference points within the reference model  
2. **Visualization:** Creates 2D slice views (axial/sagittal/coronal) in 3D reference model  

---

## How to Use
1. **Code structure:** The Intersection and Distribution Classes can be integrated into a corresponding main script. The key parts of the code are extracted into individual functions for reuse.  
2. **Import:** Import your STL mesh into MATLAB and add the *inpolyhedron* function to your path  
3. **Run Intersection Analysis:** Compute intersections for any intersection combinations  
4. **Distribution:** In relation to reference model  
5. **Visualization:** Distributions (2D/3D Views) in reference model  

---

## Use Cases
- **3D Mesh Intersection:** Determine intersection regions of multiple triangulated surfaces or volumes  
- **Frequency analysis:** Identify points that appear in multiple intersections  
- **Visualization:** Explore mesh intersections in 2D slice views (axial/sagittal/coronal)  

---

## Benefits
- **Flexible:** Works with any STL-based geometry  
- **Efficient:** Uses index masks instead of large logical arrays  
- **Scalable Intersections:** With pairwise intersections and the linked frequency tables, any intersection can be derived  
- **Evaluation:** Supports multiple intersection types (alpha, refinedAlpha, inside, grid) and compares them  

---

## Keywords
`STL`, `mesh processing`, `intersection analysis`, `overlap analysis`, `point-in-polyhedron`, `inpolyhedron`, `slice view`, `pelvis`

---

## Compatibility
These classes are compatible with MATLAB and can easily be integrated into existing workflows for analysing the intersections of anatomical models.  
It is optimised for **MATLAB 2024a** but should be compatible with the most recent versions of MATLAB.

---

## Disclaimer
This MATLAB class is provided on the MATLAB File Exchange for educational and research purposes.  
Users should ensure that the class meets their specific analysis requirements and may need to adapt it accordingly.  
The code is provided *"as-is,"* and the author assumes no responsibility for its use or any consequences thereof.

---

## References
- Johannes Korsawe (2025). [Minimal Bounding Box](https://www.mathworks.com/matlabcentral/fileexchange/18264-minimal-bounding-box), MATLAB Central File Exchange.  
- Sven (2025). [inpolyhedron - are points inside a triangulated volume?](https://www.mathworks.com/matlabcentral/fileexchange/37856-inpolyhedron-are-points-inside-a-triangulated-volume), MATLAB Central File Exchange.  

---

## Zitieren als
Carina Micheler (2025). [Anatomy IntersectionDistribution Analysis](https://www.mathworks.com/matlabcentral/fileexchange/181297-anatomy-intersectiondistribution-analysis), MATLAB Central File Exchange. Abgerufen 1. September 2025.
