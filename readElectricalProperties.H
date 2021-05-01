IOdictionary constitutiveDict
(
    IOobject
    (
        "constitutiveProperties",    // dictionary name
        runTime.constant(),     // dict is found in "constant"
        mesh,                   // registry for the dict
        IOobject::MUST_READ,    // must exist, otherwise failure
        IOobject::NO_WRITE      // dict is only read by the solver
    )
);

dimensionedScalar eps1("eps", dimensionSet(-1, -3, 4, 0, 0, 2, 0), scalar(0.0));
dimensionedScalar eps2("eps", dimensionSet(-1, -3, 4, 0, 0, 2, 0), scalar(0.0));

dimensionedScalar sig1("sig", dimensionSet(-1, -3, 3, 0, 0, 2, 0), scalar(0.0));
dimensionedScalar sig2("sig", dimensionSet(-1, -3, 3, 0, 0, 2, 0), scalar(0.0));


wordList phaseNames(constitutiveDict.lookup("phases"));

dictionary& subDict = constitutiveDict.subDict(phaseNames[0]).subDict("parameters");

eps1_ = subDict.lookup("eps");
sig1_ = subDict.lookup("sig");

dictionary& subDict = constitutiveDict.subDict(phaseNames[1]).subDict("parameters");

eps2_ = subDict.lookup("eps");
sig2_ = subDict.lookup("sig");


