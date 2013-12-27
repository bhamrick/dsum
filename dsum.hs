module Main where

import Data.Complex

type Cpx = Complex Double
type ByteDist = [Double] -- Assumed to have size 256
type ByteFFT = [Cpx] -- Fourier basis, assumed to IDFT to an all real vector

-- Encounter chances, slots are 1 indexed for simplicity
encounterThresholds = [0, 51, 102, 141, 166, 191, 216, 229, 242, 253, 256]
uniformEncounterSlot k = normalize . fromFunction $ (\x -> if x >= encounterThresholds !! (k-1) &&
                                                              x < encounterThresholds !! k
                                                           then 1 else 0)
-- Given a distribution over RNG2, compute the distribution over encounters (again 1 indexed)
encounterProbabilities :: ByteDist -> [Double]
encounterProbabilities d = threshSums encounterThresholds (zip [0..] d)
    where threshSums _ [] = []
          threshSums (t:ts) v = sum (map snd (takeWhile (\(i, p) -> i < t) v)) : threshSums ts (dropWhile (\(i, p) -> i < t) v)

framesPerStep = 17
encounterRate = 15

encounterDsumDelta = gaussian (256 - 25.39) 10.176
frameDsumDelta = fromFunction (\x -> case x of 255 -> 0.086
                                               0   -> 0.501
                                               1   -> 0.413
                                               _   -> 0)

goalSlots = [3]
successProbability :: [Int] -> ByteDist -> Double
successProbability goal d = sum (map (encounterProbabilities d !!) goal)

chancesFromSlot :: Int -> [Double]
chancesFromSlot k = map (successProbability goalSlots . (negateDist encounterChance !+ encounterChance !+ uniformEncounterSlot k !+ encounterDsumDelta !+)) stepMixers
    where stepMixers = map (\n -> frameDsumDelta !* (16 * n)) [0..]
          encounterChance = uniformRange 0 (encounterRate - 1)

main :: IO ()
main = do
    mapM_ handleSlot [1..10]
    
handleSlot :: Int -> IO ()
handleSlot k = do
    print k
    printStrategy . strategy . take 250 . chancesFromSlot $ k
    putStrLn ""

strategy :: [Double] -> [Bool]
strategy = map (>= 0.3)

printStrategy :: [Bool] -> IO ()
printStrategy = printStrategy' False where
    printStrategy' _ [] = return ()
    printStrategy' b s = do
        putStr $ (show . length . takeWhile (== b) $ s) ++ " "
        printStrategy' (not b) (dropWhile (== b) s)

-- Below here are the backend things
sep :: [a] -> ([a], [a])
sep [] = ([], [])
sep (x:[]) = ([x], [])
sep (x1:x2:xs) = let (l1, l2) = sep xs in (x1:l1, x2:l2)

rootOfUnity :: Int -> Double -> Cpx
rootOfUnity n k = cis(2*pi*k/(fromIntegral n))

-- Only works for length being a power of 2
fft :: [Cpx] -> [Cpx]
fft [] = []
fft v@(x:[]) = v
fft v = let (evens, odds) = sep v
            n = length v -- This is probably a lot of extra work
            yeven = fft evens
            yodd = fft odds
            half1 = map (\(j, ye, yo) -> ye + rootOfUnity n j * yo) (zip3 [0..] yeven yodd)
            half2 = map (\(j, ye, yo) -> ye - rootOfUnity n j * yo) (zip3 [0..] yeven yodd)
            in half1 ++ half2

fft' :: [Cpx] -> [Cpx]
fft' [] = []
fft' v@(x:[]) = v
fft' v = let (evens, odds) = sep v
             n = length v -- This is probably a lot of extra work
             yeven = fft' evens
             yodd = fft' odds
             half1 = map (\(j, ye, yo) -> ye + rootOfUnity n (-j) * yo) (zip3 [0..] yeven yodd)
             half2 = map (\(j, ye, yo) -> ye - rootOfUnity n (-j) * yo) (zip3 [0..] yeven yodd)
             in half1 ++ half2

ifft :: [Cpx] -> [Cpx]
ifft v = map (/ fromIntegral (length v)) (fft' v)

(.*) :: Num a => [a] -> [a] -> [a]
(.*) = (map (uncurry (*)) .) . zip

(.^) :: (Num a, Integral b) => [a] -> b -> [a]
v .^ n = map (^ n) v

convolve v1 v2 = ifft (fft v1 .* fft v2)

realify :: [Cpx] -> [Double]
realify = map magnitude

complexify :: [Double] -> [Cpx]
complexify = map (:+ 0)

-- Adding distributions
infixl 6 !+
(!+) :: ByteDist -> ByteDist -> ByteDist
v1 !+ v2 = realify (complexify v1 `convolve` complexify v2)

-- Add a distribution to itself multiple times
infixl 7 !*
(!*) :: Integral b => ByteDist -> b -> ByteDist
v !* b = realify . ifft $ fft (complexify v) .^ b

-- Basic distributions
fromFunction :: (Int -> Double) -> ByteDist
fromFunction f = map f [0..255]

uniform :: ByteDist
uniform = fromFunction (const $ 1/256)

uniformRange :: Int -> Int -> ByteDist
uniformRange lo hi = if lo <= hi
    then let size = fromIntegral (hi - lo + 1)
             f = (\x -> if x >= lo && x <= hi
                        then 1/size
                        else 0) in
         fromFunction f
    else let size = fromIntegral (256 - lo + hi + 1)
             f = (\x -> if x >= lo || x <= hi
                        then 1/size
                        else 0) in
             fromFunction f

binary :: Double -> ByteDist
binary p = fromFunction (\x -> case x of 0 -> (1 - p)
                                         1 -> p
                                         _ -> 0)

negateDist :: ByteDist -> ByteDist
negateDist v = head v : reverse (tail v)

normalize :: [Double] -> ByteDist
normalize v = map (/ sum v) v

gaussian :: Double -> Double -> ByteDist
gaussian mu sigma = normalize . fromFunction $ \x -> exp (-(arcDiff (fromIntegral x) mu)^2 / (2 * sigma^2))
    where arcDiff a b = min (abs (b-a)) (256 - abs (b - a))
