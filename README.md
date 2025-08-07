Polygon-Based Graph Visualization (GPU Version)

This project provides a GPU-accelerated method to visualize graphs using polygonal regions. Below is a detailed guide for installation and execution.

🔧 Installation

Note: This project must be run on Linux and requires CUDA to be installed on your system.

Clone the Repository

\`\`\`sh

git clone https://github.com/vga-usask/Map-Visualizations-For-Graphs.git 

cd PolygonBasedGraphVisualization

\`\`\`

Install CUDA (if not already installed)

Make sure CUDA is installed and configured correctly for your GPU. You can check with:

\`\`\`sh

nvcc \--version

\`\`\`

Follow the official CUDA installation guide if needed.

Build the GPU-based Code

Navigate to the build directory and compile the project:

\`\`\`sh

cd builds/linux

make graph\_viewer

\`\`\`

▶️ Running the Program

After successful compilation, you can run the program using a command like:

\`\`\`sh

./graph\_viewer ./paper/our3/small\_graph\_size\_10\_graph\_1/outputs\_convert/Vpolygon\_points.txt ./paper/our3/small\_graph\_size\_10\_graph\_1/outputs\_convert/edgelist.txt ./paper/our3/small\_graph\_size\_10\_graph\_1 2000 1 wg 1 1 approximate 1 1 1 png svg gpu 1866 ./paper/our3/small\_graph\_size\_10\_graph\_1/outputs\_convert/node\_polygon\_mapping.txt .7 1 1 ./polygon/default\_node\_map\_dep\_map.txt ./polygon/default\_node\_map\_dep\_map.txt

\`\`\`

📦 Breakdown of Command-Line Arguments

Argument

Description

argv\[1\]

Path to the polygon vertices file

argv\[2\]

Path to the edge list file

argv\[3\]

Output directory

argv\[4\]

Maximum number of iterations for the layout algorithm

argv\[5\]

Number of screenshots to generate

argv\[6\]

Gravity mode: wg (weak gravity) or sg (strong gravity)

argv\[7\]

Repulsion force strength (scale)

argv\[8\]

Gravity multiplier (typically set to 1\)

argv\[9\]

Repulsion method: approximate or exact

argv\[10\]

Repulsion force tuning flag (1 to enable)

argv\[11\]

Corner gravity activation (1 \= on)

argv\[12\]

Polygon size scaling factor

argv\[13\]

Output as PNG? (png to enable)

argv\[14\]

Output as SVG? (svg to enable)

argv\[15\]

Use GPU? (gpu to enable)

argv\[16\]

Total number of nodes in the input graph

argv\[17\]

Path to the node-to-polygon mapping file

argv\[18\]

Maximum allowed area difference for layout tuning

argv\[19\]

Exceptional polygon ID (for specific scaling treatment)

argv\[20\]

Scale adjustment for the exceptional polygon

📂 Input Files

Vpolygon\_points.txt

Each line represents a polygon using the format:

\<polygon\_id\> x1,y1 x2,y2 x3,y3 ... xn,yn

Example:

0 0.0,0.0 1.0,0.0 1.0,1.0 0.0,1.0 1 2.0,2.0 3.0,2.0 3.0,3.0 2.0,3.0

edgelist.txt

Contains the list of edges, one per line:

\<source\_node\_id\> \<target\_node\_id\>

Example:

0 1 1 2 2 3 3 0

node\_polygon\_mapping.txt

Maps node IDs to polygon IDs:

\<polygon\_id\> \<node\_id\_1\> \<node\_id\_2\> ... \<node\_id\_n\>

Example:

0 0 1 1 2 3

🖼️ Output

The program generates layout visualizations in the specified formats (png, svg) and saves them in a timestamped directory within the output path.

🔍 Reproducing Figure 1 (GI 2025 Paper)

To recreate Figure 1 from our paper "Map Visualizations for Graphs with Group Restrictions" presented at GI 2025, follow these steps:

Unzip “builds.zip”

Use the following files found in paper/our2/dlpb\_small/outputs\_convert:

polyline\_points.txt

edgelist.txt – network edge list

node\_polygon\_mapping.txt – node-to-polygon mapping

The full path to these files can be used, or move them into builds/linux.

Navigate to builds/linux and run the following command:

./graph\_viewer ./paper/our2/dlpb\_small/outputs\_convert/polyline\_point.txt ./paper/our2/dlpb\_small/outputs\_convert/edgelist.txt ./paper/our2/dlpb\_small 2000 1 wg 1 1 approximate 1 1 1 png svg gpu 26325 ./paper/our2/dlpb\_small/outputs\_convert/node\_polygon\_mapping.txt .7 1 1 ./polygon/default\_node\_map\_dep\_map.txt ./polygon/default\_node\_map\_dep\_map.txt 

View the Result:

The output will be saved as an SVG and PNG under ./output/.

Example results are shown below:

📌 Notes

The layout engine uses GPU-accelerated ForceAtlas2, tuned for subgraph fitting within polygon boundaries.

The program dynamically adjusts scaling to minimize area difference between graph convex hulls and polygon shapes.

Corner gravity helps stretch layouts to better fit polygon shapes, improving aesthetic results.  
