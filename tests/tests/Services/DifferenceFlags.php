<?php

namespace Unitmp\TmdetTest\Services;

class DifferenceFlags {
    public bool $tmpAttributes = false;
    public bool $numberOfTmRegions = false;
    public bool $newXmlHasLoop = false;
    public bool $oldXmlHasLoop = false;
    public bool $regionBoundaries = false;
    public bool $regionStart = false;
    public bool $regionEnd = false;
    // public bool $regionTypes = false;
    public bool $chainLists = false;
    public bool $addedChainLists = false;
    public bool $deletedChainLists = false;
    public bool $onlyNewTmdetXmlHasAddedChain = false;
    public bool $onlyOldTmdetXmlHasAddedChain = false;
    public bool $onlyNewTmdetXmlHasDeletedChain = false;
    public bool $onlyOldTmdetXmlHasDeletedChain = false;
    public bool $fatalErrorOccured = false;
    public bool $newRunWithNoForceDelete = false;
}
