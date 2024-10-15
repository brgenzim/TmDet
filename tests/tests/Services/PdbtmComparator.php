<?php

namespace Unitmp\TmdetTest\Services;

class PdbtmComparator {

    const MAX_ALLOWED_RESIDUE_DIFFERENCE = 4;

    public static function compareTmdetData(Tmdet30Runner $runner): array {

        $details = [];
        $messages = [];
        $categories = [];
        $flags = new DifferenceFlags();
        $flags->newRunWithNoForceDelete = $runner->noForceDelete;

        $oldTmp = $runner->oldData['isTmp'];
        $newTmp = $runner->newData['isTmp'];
        if ($oldTmp != $newTmp) {
            $category = 'TMP attributes differ';
            $flags->tmpAttributes = true;
            return [ 'categories' => [ $category ], 'flags' => $flags ];
        }

        $oldDeleteChains = $runner->oldData['deletedChains'];
        $newDeleteChains = $runner->newData['deletedChains'];
        $oldApplyToChains = $runner->oldData['addedChains'];
        $newApplyToChains = $runner->newData['addedChains'];

        static::updateDetails(
            condition: empty($newDeleteChains) && !empty($oldDeleteChains), flags: $flags, flagName: 'onlyOldTmdetXmlHasDeletedChain',
            categories: $categories, category: 'Old TMDET XML has deleted chains, but new XML hasn\'t contain any'
        );

        static::updateDetails(
            condition: !empty($newDeleteChains) && empty($oldDeleteChains), flags: $flags, flagName: 'onlyNewTmdetXmlHasDeletedChain',
            categories: $categories, category: 'New TMDET XML has deleted chains, but old XML hasn\'t contain any'
        );

        static::updateDetails(
            condition: $oldDeleteChains != $newDeleteChains, flags: $flags, flagName: 'deletedChainLists',
            categories: $categories, category: 'Deleted chain list differs'
        );

        static::updateDetails(
            condition: empty($newApplyToChains) && !empty($oldApplyToChains), flags: $flags, flagName: 'onlyOldTmdetXmlHasAddedChain',
            categories: $categories, category: 'Old TMDET XML has added chains, but new XML hasn\'t contain any'
        );
        static::updateDetails(
            condition: !empty($newApplyToChains) && empty($oldApplyToChains), flags: $flags, flagName: 'onlyNewTmdetXmlHasAddedChain',
            categories: $categories, category: 'New TMDET XML has added chains, but old XML hasn\'t contain any'
        );

        static::updateDetails(
            condition: $oldApplyToChains != $newApplyToChains, flags: $flags, flagName: 'addedChainLists',
            categories: $categories, category: 'Added chain list differs'
        );

        $oldChains = $runner->oldData['chains'];
        $newChains = $runner->newData['chains'];
        $details['oldData'] = [
            'isTmp' => $oldTmp,
            'chains' => array_keys($oldChains),
            'addedChains' => $oldApplyToChains,
            'deletedChains' => $oldDeleteChains
        ];
        $details['newData'] = [
            'isTmp' => $newTmp,
            'chains' => array_keys($newChains),
            'addedChains' => $newApplyToChains,
            'deletedChains' => $newDeleteChains
        ];
        // Undo biomatrix operations and/or synchronize the different chain sets
        static::synchronizeBioMatrixDeletes($oldChains, $newChains, $oldDeleteChains, $newDeleteChains);
        static::synchronizeBioMatrixApply($oldChains, $oldApplyToChains);

        static::updateDetails(
            condition: !empty(array_diff_key($oldChains, $newChains)) || !empty(array_diff_key($newChains, $oldChains)),
            flags: $flags, flagName: 'chainLists',
            categories: $categories, category: 'Different chain lists in TMDET results'
        );

        $helixFilter = function ($regions) {
            return array(...array_filter($regions, function ($value) {
                return $value['type'] === 'H';
            }));
        };
        $loopFilter = function ($regions) {
            return array(...array_filter($regions, function ($value) {
                return $value['type'] === 'L';
            }));
        };
        $newChainHasLoop = false;
        $oldChainHasLoop = false;
        foreach ($oldChains as $chain => $chainData) {
            $category = 'Number of TM regions differs';
            $message = $category . ' - chain: ' . $chain;
            static::updateDetails(
                condition: $chainData['num_tm'] != $newChains[$chain]['num_tm'],
                flags: $flags, flagName: 'numberOfTmRegions',
                messages: $messages, categories: $categories,
                message: $message, category: $category
            );
            $newRegions = $newChains[$chain]['regions'];
            $oldRegions = $chainData['regions'];
            if (!empty($loopFilter($newRegions))) {
                $newChainHasLoop = true;
            }
            if (!empty($loopFilter($oldRegions))) {
                $oldChainHasLoop = true;
            }
            $category = 'New XML has loop';
            $message = $category . ' - chain: ' . $chain;
            static::updateDetails(
                condition: $newChainHasLoop,
                // NOTE: related flags will be set after this 'foreach'
                messages: $messages, categories: $categories,
                message: $message, category: $category
            );
            $category = 'Old XML has loop';
            $message = $category . ' - chain: ' . $chain;
            static::updateDetails(
                condition: $oldChainHasLoop,
                // NOTE: related flags will be set after this 'foreach'
                messages: $messages, categories: $categories,
                message: $message, category: $category
            );

            //
            // Compare region boundaries of helicies
            $newRegions = $helixFilter($newRegions);
            foreach ($helixFilter($oldRegions) as $index => $oldRegion) {
                // check helices are at the same index
                if (!array_key_exists($index, $newRegions)) {
                    continue;
                }
                $newRegion = $newRegions[$index];

                // NOTE: Inactivated since only helices verified from old and new regions
                // if ($oldRegion['type'] == 'H') {
                //     static::updateDetails(
                //         condition: $oldRegion['type'] != $newRegion['type'],
                //         flags: $flags, flagName: 'regionTypes',
                //         messages: $messages, categories: $categories,
                //         message: 'Region type of chain "' . $chain . '"[new seq_beg: ' . $newRegion['seq_beg'] . '] differs',
                //         category: 'Region type of chain differs'
                //     );
                // }

                $actualValue = $newRegion['seq_beg'];
                $condition = !static::equalsWithDelta($oldRegion['seq_beg'], $actualValue, static::MAX_ALLOWED_RESIDUE_DIFFERENCE);
                static::updateDetails(
                    condition: $condition,
                    messages: $messages, categories: $categories,
                    message: 'Region start in chain \'' . $chain . '\'[new seq_beg: ' . $actualValue . '] differs',
                    category: 'Region start or end in chain differs'
                );
                $flags->regionBoundaries |= $condition;

                $actualValue = $newRegion['seq_end'];
                $condition = !static::equalsWithDelta($oldRegion['seq_end'], $actualValue, static::MAX_ALLOWED_RESIDUE_DIFFERENCE);
                static::updateDetails(
                    condition: $condition,
                    messages: $messages, categories: $categories,
                    message: 'Region end in chain \'' . $chain . '\'[new seq_end: ' . $actualValue . '] differs',
                    category: 'Region start or end in chain differs'
                );
                $flags->regionBoundaries |= $condition;
            }
        }
        $flags->newXmlHasLoop = $newChainHasLoop;
        $flags->oldXmlHasLoop = $oldChainHasLoop;
        $flags->regionStart = !empty(preg_grep('/Region start in chain/', $messages));
        $flags->regionEnd = !empty(preg_grep('/Region end in chain/', $messages));

        $details['flags'] = $flags;
        $details['categories'] = $categories;
        if (!empty($messages)) {
            $details['messages'] = $messages;
        }
        return $details;
    }

    protected static function updateDetails(bool $condition, array &$categories, string $category,
        DifferenceFlags|false $flags = false, string $flagName = 'N/A',
        array|false &$messages = false, string|false $message = false): void {

        if ($flags) {
            $flags->{$flagName} = $condition;
        }

        static::addToCategory($condition, $categories, $category);
        if (!$message || $messages === false) {
            return;
        }
        static::appendMessage($condition, $messages, $message);
    }

    protected static function appendMessage(bool $condition, array &$messages, string $message): void {
        if ($condition) {
            $messages[] = $message;
        }
    }

    protected static function addToCategory(bool $condition, array &$categories, string $category): void {
        if ($condition && array_search($category, $categories) === false) {
            $categories[] = $category;
        }
    }

    public static function equalsWithDelta(int $expected, int $actual, int $delta): bool {
        return abs($expected - $actual) <= $delta;
    }

    public static function synchronizeBioMatrixDeletes(array &$oldChains, array &$newChains, array $oldDeleteChains, array $newDeleteChains) {
        foreach (array_keys($oldChains) as $chainId) {
            if (array_search($chainId, $newDeleteChains) !== false) {
                unset($oldChains[$chainId]);
            }
        }
        foreach (array_keys($newChains) as $chainId) {
            if (array_search($chainId, $oldDeleteChains) !== false) {
                unset($newChains[$chainId]);
            }
        }
    }

    public static function synchronizeBioMatrixApply(array &$chains, array $applyChains) {
        foreach (array_keys($chains) as $chainId) {
            if (array_search($chainId, $applyChains) !== false) {
                unset($chains[$chainId]);
            }
        }
    }

}


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
