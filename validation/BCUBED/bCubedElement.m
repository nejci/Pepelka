function [F, precision, recall] = bCubedElement(refrPart, respPart)

precision = bCubedElementPrecision(refrPart, respPart);
recall = bCubedElementRecall(refrPart, respPart);
F = 2 * precision * recall / (precision + recall);