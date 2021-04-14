#ifndef COLUMN_H
#define COLUMN_H

#include "Data.h"
#include "Parameters.h"

class Column{

    public:
        Column();
        /* Column(std::vector<std::pair<int,int>> dimTagsFragmentedColumn, Parameters * prm); */
        Column(std::vector<std::pair<int,int>> dimTagsFragmentedColumn, Parameters * prm, std::string _periodic);
        virtual ~Column();

        //TODO: Consider this approach
        // Column(std::vector<std::pair<int,int>> dimTagsFragmented);
        // Column::prepare();
        // Column::write(); // to handle writing all the meshes appropriately
        // This will allow me to completely clean up the Model space, and use it as a driver for building columns.
        // So, in Model..
        //      Column(Main)
        //      Column(Inlet)
        //      Column(Outlet)
        //      LinkColumnsPeriodic(Inlet, Main)
        //      LinkColumnsPeriodic(Main, Outlet)
        //      Mesh()
        //      Write()

        std::string periodic;
        double dx, dy, dz;

        Volumes volumes;
        Surfaces surfaces;
        Walls outerWalls, beadWalls;

        void generateBoxSurfaces();
        void AddPhysicalGroups();
        void separateBoundingSurfaces(std::vector<std::pair<int, int>> dimTagsSurfaces, Walls& tWall);

        void stats();

        void setupPeriodicSurfaces(Walls& tWall);
        void matchPeriodicSurfaces(std::vector<int>& ltags, std::vector<int>& rtags, int per_dir, std::vector<double> affineTranslation);
        void linkPeriodicZ(Column& second);

        void write(std::string outpath, std::string outfile, std::string extension);
        void writeFragments(std::string outpath, std::string outfile, std::string extension);
        void meshVolumes(Parameters * prm);
};

#endif
