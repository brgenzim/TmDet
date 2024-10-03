<?php

namespace Unitmp\TmdetTest\Services;

class PdbtmComparator {

    const MAX_ALLOWED_RESIDUE_DIFFERENCE = 4;

    public static function compareTmdetData(Tmdet30Runner $runner): array {

        $details = [];
        $messages = [];

        $oldValue = $runner->oldData['isTmp'];
        $newValue = $runner->newData['isTmp'];
        if ($oldValue != $newValue) {
            static::appendMessage(true, $messages, '"isTmp" values differ');
            return [ 'messages' => $messages ];
        }

        $oldDeleteChains = $runner->oldData['deletedChains'];
        $newDeleteChains = $runner->newData['deletedChains'];
        $oldApplyToChains = $runner->oldData['addedChains'];
        $newApplyToChains = $runner->newData['addedChains'];

        $errorEncountered = empty($newDeleteChains) && !empty($oldDeleteChains);
        static::appendMessage($errorEncountered, $messages, 'Old TMDET XML has deleted chains, but new XML hasn\'t contain any');
        $errorEncountered = !empty($newDeleteChains) && empty($oldDeleteChains);
        static::appendMessage($errorEncountered, $messages, 'New TMDET XML has deleted chains, but old XML hasn\'t contain any');
        static::appendMessage($oldDeleteChains != $newDeleteChains, $messages, 'Deleted chain list differs');

        $errorEncountered = empty($newApplyToChains) && !empty($oldApplyToChains);
        static::appendMessage($errorEncountered, $messages, 'Old TMDET XML has added chains, but new XML hasn\'t contain any');
        $errorEncountered = !empty($newApplyToChains) && empty($oldApplyToChains);
        static::appendMessage($errorEncountered, $messages, 'New TMDET XML has added chains, but old XML hasn\'t contain any');
        static::appendMessage($oldApplyToChains != $newApplyToChains, $messages, 'Added chain list differs');

        $oldChains = $runner->oldData['chains'];
        $newChains = $runner->newData['chains'];
        $details['oldData'] = [
            'chains' => array_keys($oldChains),
            'addedChains' => $oldApplyToChains,
            'deletedChains' => $oldDeleteChains
        ];
        $details['newData'] = [
            'chains' => array_keys($newChains),
            'addedChains' => $newApplyToChains,
            'deletedChains' => $newDeleteChains
        ];
        // Undo biomatrix operations and/or synchronize the different chain sets
        static::synchronizeBioMatrixDeletes($oldChains, $newChains, $oldDeleteChains, $newDeleteChains);
        static::synchronizeBioMatrixApply($oldChains, $oldApplyToChains);

        $errorEncountered = !empty(array_diff_key($oldChains, $newChains)) || !empty(array_diff_key($newChains, $oldChains));
        static::appendMessage($errorEncountered, $messages, 'Different chain lists in TMDET results');

        $helixFilter = function ($regions) {
            return array_filter($regions, function ($value) {
                return $value === 'H';
            } );
        };
        foreach ($oldChains as $chain => $chainData) {
            static::appendMessage($chainData['num_tm'] != $newChains[$chain]['num_tm'], $messages,
                'Number of TM regions differs - chain: ' . $chain);
            $newRegions = $helixFilter($newChains[$chain]['regions']);
            foreach ($helixFilter($chainData['regions']) as $index => $oldRegion) {
                $newRegion = $newRegions[$index];
                if ($oldRegion['type'] == 'H') {
                    static::appendMessage($oldRegion['type'] != $newRegion['type'], $messages,
                        'Region type of chain "' . $chain . '"[new seq_beg: ' . $newRegion['seq_beg'] . '] differs');
                }

                $actualValue = $newRegion['seq_beg'];
                $errorEncountered = !static::equalsWithDelta($oldRegion['seq_beg'], $actualValue, static::MAX_ALLOWED_RESIDUE_DIFFERENCE);
                static::appendMessage($errorEncountered, $messages,
                    'Region start in chain "' . $chain . '"[new seq_beg: ' . $actualValue . '] differs');

                $actualValue = $newRegion['seq_end'];
                $errorEncountered = !static::equalsWithDelta($oldRegion['seq_end'], $actualValue, static::MAX_ALLOWED_RESIDUE_DIFFERENCE);
                static::appendMessage($errorEncountered, $messages,
                    'Region end in chain "' . $chain . '"[new seq_end: ' . $actualValue . '] differs');
            }
        }

        $details['messages'] = $messages;
        return $details;
    }

    protected static function appendMessage(bool $condition, array &$messages, string $message): void {
        if ($condition) {
            $messages[] = $message;
        }
    }

    public static function equalsWithDelta(int $expected, int $actual, int $delta): bool {
        return abs($expected - $actual) < $delta;
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