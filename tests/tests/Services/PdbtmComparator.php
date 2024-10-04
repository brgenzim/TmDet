<?php

namespace Unitmp\TmdetTest\Services;

class PdbtmComparator {

    const MAX_ALLOWED_RESIDUE_DIFFERENCE = 4;

    public static function compareTmdetData(Tmdet30Runner $runner): array {

        $details = [];
        $messages = [];
        $categories = [];

        $oldTmp = $runner->oldData['isTmp'];
        $newTmp = $runner->newData['isTmp'];
        if ($oldTmp != $newTmp) {
            $category = 'TMP attributes differ';
            return [ 'categories' => [ $category ] ];
        }

        $oldDeleteChains = $runner->oldData['deletedChains'];
        $newDeleteChains = $runner->newData['deletedChains'];
        $oldApplyToChains = $runner->oldData['addedChains'];
        $newApplyToChains = $runner->newData['addedChains'];

        static::updateDetails(
            condition: empty($newDeleteChains) && !empty($oldDeleteChains),
            categories: $categories, category: 'Old TMDET XML has deleted chains, but new XML hasn\'t contain any'
        );

        static::updateDetails(
            condition: !empty($newDeleteChains) && empty($oldDeleteChains),
            categories: $categories, category: 'New TMDET XML has deleted chains, but old XML hasn\'t contain any'
        );

        static::updateDetails(
            condition: $oldDeleteChains != $newDeleteChains,
            categories: $categories, category: 'Deleted chain list differs'
        );

        static::updateDetails(
            condition: empty($newApplyToChains) && !empty($oldApplyToChains),
            categories: $categories, category: 'Old TMDET XML has added chains, but new XML hasn\'t contain any'
        );
        static::updateDetails(
            condition: !empty($newApplyToChains) && empty($oldApplyToChains),
            categories: $categories, category: 'New TMDET XML has added chains, but old XML hasn\'t contain any'
        );

        static::updateDetails(
            condition: $oldApplyToChains != $newApplyToChains,
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
            categories: $categories, category: 'Different chain lists in TMDET results'
        );

        $helixFilter = function ($regions) {
            return array_filter($regions, function ($value) {
                return $value === 'H';
            } );
        };
        foreach ($oldChains as $chain => $chainData) {
            $category = 'Number of TM regions differs';
            $message = $category . ' - chain: ' . $chain;
            static::updateDetails(
                condition: $chainData['num_tm'] != $newChains[$chain]['num_tm'],
                messages: $messages, categories: $categories,
                message: $message, category: $category
            );
            $newRegions = $helixFilter($newChains[$chain]['regions']);
            foreach ($helixFilter($chainData['regions']) as $index => $oldRegion) {
                $newRegion = $newRegions[$index];
                if ($oldRegion['type'] == 'H') {
                    static::updateDetails(
                        condition: $oldRegion['type'] != $newRegion['type'],
                        messages: $messages, categories: $categories,
                        message: 'Region type of chain "' . $chain . '"[new seq_beg: ' . $newRegion['seq_beg'] . '] differs',
                        category: 'Region type of chain differs'
                    );
                }

                $actualValue = $newRegion['seq_beg'];
                static::updateDetails(
                    condition: !static::equalsWithDelta($oldRegion['seq_beg'], $actualValue, static::MAX_ALLOWED_RESIDUE_DIFFERENCE),
                    messages: $messages, categories: $categories,
                    message: 'Region start in chain "' . $chain . '"[new seq_beg: ' . $actualValue . '] differs',
                    category: 'Region start/end in chain differs'
                );

                $actualValue = $newRegion['seq_end'];
                static::updateDetails(
                    condition: !static::equalsWithDelta($oldRegion['seq_end'], $actualValue, static::MAX_ALLOWED_RESIDUE_DIFFERENCE),
                    messages: $messages, categories: $categories,
                    message: 'Region end in chain "' . $chain . '"[new seq_end: ' . $actualValue . '] differs',
                    category: 'Region start/end in chain differs'
                );
            }
        }

        $details['categories'] = $categories;
        if (!empty($messages)) {
            $details['messages'] = $messages;
        }
        return $details;
    }

    protected static function updateDetails(bool $condition, array &$categories, string $category,
        array|false &$messages = false, string|false $message = false): void {

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