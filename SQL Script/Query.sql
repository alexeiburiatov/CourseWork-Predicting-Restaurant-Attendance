CREATE DATABASE RestaurantVisitingDB
USE RestaurantVisitingDB

DROP TABLE IF EXISTS Dim_Cuisine
CREATE TABLE Dim_Cuisine(
cuisineID INT PRIMARY KEY NOT NULL,
cuisineName NVARCHAR(255 ) NOT NULL,
loadDate DATETIME DEFAULT(GETDATE()) NOT NULL,
loadEndDate DATETIME DEFAULT('9999-01-01 23:59:59') NOT NULL
)

DROP TABLE IF EXISTS Dim_Place
CREATE TABLE Dim_Place(
placeID INT PRIMARY KEY NOT NULL,
restaurantName NVARCHAR(255) NOT NULL,
cuisineID INT NOT NULL,
latitude NVARCHAR(255) NOT NULL,
longitude NVARCHAR(255) NOT NULL,
CONSTRAINT FK_Dim_Cuisine_Dim_Place FOREIGN KEY (cuisineID)
    REFERENCES Dim_Cuisine (cuisineID)
    ON DELETE CASCADE
	ON UPDATE CASCADE
)


DROP TABLE IF EXISTS Dim_Reservation
CREATE TABLE Dim_Reservation(
reserveID INT PRIMARY KEY NOT NULL,
dateID INT NOT NULL,
timeValue TIME NOT NULL,
resrvAmntOfVisitors INT NOT NULL,
placeID INT NOT NULL,
CONSTRAINT FK_Dim_Date_Dim_Reservation FOREIGN KEY (dateID)
    REFERENCES Dim_Date (dateID)
    ON DELETE CASCADE
	ON UPDATE CASCADE,
CONSTRAINT FK_Dim_Place_Dim_Reservation FOREIGN KEY (placeID)
    REFERENCES Dim_Place (placeID)
    ON DELETE CASCADE
	ON UPDATE CASCADE
)

CREATE TABLE Dim_DayOfTheWeek(
dayOfTheWeekID INT PRIMARY KEY NOT NULL,
dayOfTheWeekName NVARCHAR(255) NOT NULL
)

DROP TABLE IF EXISTS Dim_Date
CREATE TABLE Dim_Date(
dateID INT PRIMARY KEY NOT NULL,
dayValue DATE NOT NULL,
dayOfTheWeekID INT NOT NULL,
isHoliday BIT NOT NULL,
CONSTRAINT FK_Dim_DayOfTheWeek_Dim_Date FOREIGN KEY (dayOfTheWeekID)
    REFERENCES Dim_DayOfTheWeek (dayOfTheWeekID)
    ON DELETE CASCADE
	ON UPDATE CASCADE

)

DROP TABLE IF EXISTS Fact_RestaurantVisiting
CREATE TABLE Fact_RestaurantVisiting(
dateID INT NOT NULL,
placeID INT NOT NULL,
amountOfVisitors INT NOT NULL,
CONSTRAINT FK_Dim_Date_Fact_RestaurantVisiting FOREIGN KEY (dateID)
    REFERENCES Dim_Date (dateID)
    ON DELETE CASCADE
	ON UPDATE CASCADE,
CONSTRAINT FK_Dim_Place_Fact_RestaurantVisiting FOREIGN KEY (placeID)
    REFERENCES Dim_Place (placeID)
    ON DELETE CASCADE
	ON UPDATE CASCADE,
CONSTRAINT PK_dateID_placeID PRIMARY KEY (dateID, placeID)

CREATE OR ALTER FUNCTION dbo.getFactAmountVisitors(
@placeName NVARCHAR(255)
)
RETURNS TABLE
AS RETURN (
	SELECT dd.dayValue
		,frv.amountOfVisitors
		,wd.dayOfTheWeekName
		,dd.isHoliday
	FROM dbo.Fact_RestaurantVisiting frv
	JOIN dbo.Dim_Date dd ON dd.dateID=frv.dateID
	JOIN dbo.Dim_DayOfTheWeek wd ON wd.dayOfTheWeekID=dd.dayOfTheWeekID
	WHERE frv.placeID=(
	SELECT dp.placeID FROM dbo.Dim_Place dp
	WHERE dp.restaurantName=@placeName
	)
);
GO


CREATE OR ALTER FUNCTION dbo.getFactAmountResrvVisitors(
@placeName NVARCHAR(255)
)
RETURNS TABLE
AS RETURN(
	SELECT DISTINCT dd.dayValue
		,SUM(dr.resrvAmntOfVisitors) OVER(PARTITION BY dd.dayValue) AS reservSum
		,wd.dayOfTheWeekName
		,dd.isHoliday
	FROM dbo.Dim_Reservation dr
	JOIN dbo.Dim_Date dd ON dd.dateID=dr.dateID
	JOIN dbo.Dim_DayOfTheWeek wd ON wd.dayOfTheWeekID=dd.dayOfTheWeekID
	WHERE dr.placeID =(
	SELECT dp.placeID FROM dbo.Dim_Place dp
	WHERE dp.restaurantName=@placeName
	)
	);
GO

SELECT * FROM dbo.getFactAmountResrvVisitors('Hyogo-ken Kakogawa Kakogawacho Kitazaike')
