{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE TypeFamilies #-}
module Main where

import Test.Tasty
import           Test.Tasty.HUnit
import Bio.ChromVAR
import Data.List
import qualified Eigen.Matrix as E
import Bio.Data.Bed.Types
import Data.Maybe
import Conduit
import qualified Eigen.SparseMatrix as ES
import qualified Eigen.Arithmetic as A

main :: IO ()
main = defaultMain $ testGroup "Main"
    [ testCase "sortBed" test ]

test :: Assertion
test = res @=? (deviations, z)
  where
    [res] = runIdentity $ runConduit $ yieldMany [(4, cellByPeak)] .|
        computeDeviation 3 3 expectation bg peakBymotif .| sinkList
    expectation = let e = E.ones `A.mul` (ES.fromList cellByPeak :: ES.SparseMatrix 4 3 Double) :: E.Matrix 1 3 Double
                      s = E.sum e
                  in concat $ E.toList $ E.map (/s) e
    round' x = fromIntegral (round (x * 1e7)) / 1e7

cellByPeak :: [(Int, Int, Double)]
cellByPeak = [(0,1,1), (0,2,2), (1,0,1), (1,1,1), (2,0,2), (2,2,1), (3,1,1), (3,2,1)]

peakBymotif :: [(Int, Int, Double)]
peakBymotif = [(0,0,1), (0,1,1), (1,0,1), (1,1,1), (1,2,1), (2,0,1), (2,2,1)]

peaks :: [BED3]
peaks = map (\(chr, i) -> BED3 chr i (i + 500)) dat
  where
    dat = [ ("chr1", 76585873)
          , ("chr2", 42772928)
          , ("chr2", 100183786) ]

bg :: [[Int]]
bg = [ [1,0,1]
     , [1,0,1]
     , [1,2,0] ]

deviations = transpose [ [0.1728395, -0.4444444, 0.1728395, -0.07407407]
             , [-0.2910053, 0.3174603, 0.2116402, -0.19841270]
             , [0.7407407, -0.6349206, -0.7407407,  0.63492063] ]

z = transpose [ [1.1547005, -1.1547005, 1.1547005, -1.1547005]
    , [-0.5773503, 0.5773503, 0.5773503, -0.5773503]
    , [3.2331615, -1.1547005, -4.0414519, 9.2376043] ]
