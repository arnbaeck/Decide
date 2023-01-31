# Decide

## The Decide program
**Decide** is a program which decides whether to launch an interceptor based on input data and *Launch Interceptor Conditions* (LIC's). Depending on the conditions, it should print "YES" for launch, or "NO" otherwise.  

## Arrays and Vectors used
*Conditions Met Vector (CMV)* - A boolean, 15-element vector whose values indicate whether the conditions of the corresponding LIC is met.
*Logical Connector Matrix (LCM)* - A 15x15 matrix which defines which LIC's are to be considered jointly, with values ANDD if both LIC's must be true, ORR if at least one of them or NOTUSED if they are not used.
*Preliminary Unlocking Matrix (PUM)* - Is a boolean 15x15 matrix which holds the resulting combined values based on the CMV and LCM.
*Preliminary Unlocking Vector (PUV)* - Is a boolean, 15-element vector that holds which LIC's that are actually relevant in the current launch determination.
*Final Unlocking Vector (FUV)* - A boolean, 15-element vector which holds the final verdict. All 15 values have to be true to unlock the launch signal.

## Contributions

# 2023-01-27
At the first meeting, Anton, Carin and Oguz contributed in creating issues for each of the LIC's, an issue for the code skeleton and and an issue for the complete unit test. (Gabriel was unable to attend at the creation of issues, but were given issues to create in the later stages as compensation). A code skeleton was created, using live share where everyone contributed in the same file. This was to make sure everyone got familiar with the base of the code. Everyone experimented and researched with gridle to create a build and version control through GitHub. All of the LIC's where then distributed among the four members (using a program on internet to randomly distribute them).

Anton: LIC 5, 6, 9, 13. 
Carin: LIC 2, 4, 10, 12.
Oguz: 1, 3, 14.
Gabriel: 0, 7, 8, 11.

Pull requests were made for each of the method implementations, it was made sure that everyone received a roughly equal amount of issues to merge. Some additional issues were created along the way to fix bugs and add missing documentation. After Gabriel had created a few issues as compensation for missing out at the beginning, the creation, resolving and merging of these issues and implementations were also attempted to evenly distribute among the members. However, while the issues were distributed among the members, some collaboration occurred when discussing how to write or fix problems or implementations. This was done over slack during the weekend, or in person when we sat on campus on monday. 

For P+:
The most commits are linked to an issue which describes the feature.
 