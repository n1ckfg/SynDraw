import pymeshlab as ml
import latk
from pathlib import Path
import os
import xml.etree.ElementTree as etree
import latk

samplePercentage = 0.3
inputFormat = "ply"
inputPath = "input"

la = latk.Latk()
la.layers.append(latk.LatkLayer())
counter = 0
minStrokeLength = 2

for fileName in os.listdir(inputPath):
    if fileName.endswith(inputFormat): 
        frame = latk.LatkFrame()
        la.layers[0].frames.append(frame)

        inputUrl = os.path.join(inputPath, fileName)
        print("Loading " + inputUrl)

        ms = ml.MeshSet()
        ms.load_new_mesh(inputUrl)
        mesh = ms.current_mesh()

        newSampleNum = int(mesh.vertex_number() * samplePercentage)
        if (mesh.edge_number() == 0 and mesh.face_number() == 0):
            ms.poisson_disk_sampling(samplenum=newSampleNum, subsample=True)
        else:
            ms.poisson_disk_sampling(samplenum=newSampleNum, subsample=False)
        ms.surface_reconstruction_ball_pivoting()
        ms.vertex_attribute_transfer(sourcemesh=0, targetmesh=1)
        ms.save_current_mesh("input.ply", save_vertex_color=True)

        os.system("SynDraw -p test.template")

        print("parsing SVG")
        tree = etree.parse("out.svg")
        root = tree.getroot()
        for element in root:
            if (element.tag.endswith("polyline")):
                points = element.get("points3d").split(" ")
                lPoints = []
                for point in points:
                    point2 = point.split(",")
                    try:
                        point3 = (float(point2[0]), float(point2[2]), float(point2[1]))
                        lPoints.append(latk.LatkPoint(point3))
                    except:
                        pass
                if (len(lPoints) >= minStrokeLength):
                    la.layers[0].frames[counter].strokes.append(latk.LatkStroke(lPoints))
        print("Finished frame " + str(counter+1))
        counter += 1

la.write("output.latk")
print("Wrote latk")

