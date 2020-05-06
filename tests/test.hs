{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE PartialTypeSignatures #-}
module Main where

import Test.Tasty
import           Test.Tasty.HUnit
import Data.List
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Storable as VS
import Bio.Data.Bed.Types
import Data.Maybe
import Conduit
import qualified Data.Matrix.Static.Dense as D
import qualified Data.Matrix.Static.Generic as D
import qualified Data.Matrix.Static.Sparse as S
import Data.Matrix.Dynamic (Dynamic(..), matrix)
import Data.Matrix.Static.LinearAlgebra
import Statistics.Sample hiding (covariance)
import System.Random.MWC (create)

import Bio.ChromVAR
import Bio.ChromVAR.Utils

main :: IO ()
main = defaultMain $ testGroup "Main"
    [ testCase "sortBed" test
    , testCase "sortBed" bgTest]
    --, testCase "sortBed" peakTest ]

test :: Assertion
test = (map round' $ concat dev', map round' $ concat z') @=?
    (map round' $ concat deviations, map round' $ concat z)
  where
    [(dev', z')] = runIdentity $ runConduit $ yieldMany [(4, cellByPeak)] .|
        computeDeviation 3 3 expectation bg peakBymotif .| sinkList
    expectation = let e = D.replicate 1 @@ (S.fromTriplet (U.fromList cellByPeak) :: SparseMatrix 4 3 Double) :: Matrix 1 3 Double
                      s = VS.sum $ D.flatten e
                  in D.toList $ D.map (/s) e
    round' :: Double -> Double
    round' x = fromIntegral (round (x * 1e4)) / 1e4

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

bgTest :: Assertion
bgTest = do
    input <- readData "tests/data/coordinates.tsv"
    expected <- readData "tests/data/transformed.tsv"
    let actual = case matrix (map (\(a,b) -> [a,b]) $ U.toList input) of
            Dynamic mat@(D.Matrix _) -> U.fromList $ map (\[a,b] -> (a,b)) $
                map VS.toList $ D.toRows $ whiten Cholesky mat
    U.map (\(x,y) -> (round' x, round' y)) actual @=?
        U.map (\(x,y) -> (round' x, round' y)) expected
  where
    readData fl = U.fromList . map (f . words) . lines <$> readFile fl
      where
        f [x,y] = (read x, read y)
    round' :: Double -> Double
    round' x = fromIntegral (round (x * 1e4)) / 1e4

peakTest :: Assertion
peakTest = do
    coord <- readData "tests/data/coordinates.tsv"

    expected <- map (map (\x -> coord U.! (read x - 1)) . words) .
        lines <$> readFile "tests/data/background.tsv"
    actual <- fmap (transpose . map (map (coord U.!) . U.toList)) $
        create >>= getBackgroundPeaks 50 coord

    print $ zipWith f expected actual

    --print $ flip map input' $ \x -> meanVarianceUnb $ U.fromList $ map dist $ comb x
    return undefined
  where
    f xs ys = meanVarianceUnb $ U.fromList $ map dist $ map (\[a,b] -> (a,b)) $ sequence [xs, ys]
    readData fl = U.fromList . map (f . words) . lines <$> readFile fl
      where
        f [x,y] = (read x, read y)
    comb (x:xs) = zip (repeat x) xs ++ comb xs
    comb _ = []
    dist ((x1,y1), (x2,y2)) = sqrt $ (x1 - x2)**2 + (y1-y2)**2