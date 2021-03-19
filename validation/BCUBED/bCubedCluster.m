function [F, precision, recall] = bCubedCluster(refrPart, respPart)

precision = bCubedClusterPrecision(refrPart, respPart);
recall = bCubedClusterRecall(refrPart, respPart);
F = 2 * precision * recall / (precision + recall);