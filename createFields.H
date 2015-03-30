    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

     dimensionedScalar dim_s
     (
         transportProperties.lookup("dim_s")
     );
    dimensionedScalar dim_m
    (
        transportProperties.lookup("dim_m")
    );

    dimensionedScalar g
    (
        transportProperties.lookup("g")
    );
    dimensionedScalar Hdry
    (
        transportProperties.lookup("Hdry")
    );
    dimensionedScalar Hdry2
    (
        transportProperties.lookup("Hdry2")
    );

    dimensionedScalar Cnu
    (
        transportProperties.lookup("Cnu")
    );

    dimensionedScalar nutmax
    (
        transportProperties.lookup("nutmax")
    );


    Info<< "Reading field H\n" << endl;
    Info<< "Flow depth\n" << endl;
    volScalarField H
    (
        IOobject
        (
            "H",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    volScalarField Hclip
    (
        IOobject
        (
            "Hclip",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
	(H-Hdry)*pos(H-Hdry) + Hdry        
    );


    Info<< "Reading field S\n" << endl;
    Info<< "Bottom elevation\n" << endl;
    volScalarField S
    (
        IOobject
        (
            "S",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    Info<< "Reading field kst\n" << endl;
    Info<< "Strickler value\n" << endl;
    volScalarField kst
    (
        IOobject
        (
            "kst",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

   Info<< "Reading field HU\n" << endl;
    volVectorField HU
    (
        IOobject
        (
            "HU",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Creating field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        HU/dim_m
    );

    Info<< "Creating field phi\n" << endl;
    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(U) & mesh.Sf()
    );

    Info<< "Creating field alpha\n" << endl;
    volScalarField alpha
    (
        IOobject
        (
            "alpha",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        (g*dim_s*dim_s/dim_m)*mag(U)/pow(kst,2.0)/pow(Hclip/dim_m,1.0/3.0)/Hclip
    );


   Info<< "Creating field nut\n" << endl;
    volScalarField nut
    (
        IOobject
        (
            "nut",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
	);


   Info<< "Creating field wetDry\n" << endl;
    volScalarField wetDry
    (
        IOobject
        (
            "wetDry",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        pos(H-Hdry)
	);


    Info<< "createFields done\n" << endl;
