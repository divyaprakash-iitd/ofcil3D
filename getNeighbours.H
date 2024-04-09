label cellIndex = mesh.findCell(rr);
// std::cout << "Cell index: " << cellIndex << std::endl;
// 
// Info << "F: " << F[cellIndex] << endl;

// Try to make a list of list so that neighbours don't have to be searched two times.
// Once for calculating forces and once more for updating positions
// labelList appended;
DynamicList<Foam::label> appended;
appended.append(cellIndex);

// Info << "Cell Index: " << cellIndex << endl;

// Get the neighbouring cells 4 times
for (int l = 0; l < 4; ++l) {
    int nn = appended.size();
    for (int m = 0; m < nn; ++m) {
        labelList neighbors = mesh.cellCells()[appended[m]];
        // appended.append(neighbors);
        appended.append(neighbors);
    }
    // Get rid of repeating values and also sort the array
    Foam::inplaceUniqueSort(appended);
}


// Info << "List: " << appended << endl;
// Info << "ListofLists: " << neighborsList[0] << endl;
