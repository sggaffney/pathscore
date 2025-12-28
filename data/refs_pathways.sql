-- PathScore Pathways Table
-- This is a minimal schema for the pathways table in the refs database.
-- Real pathway data should be imported separately.

CREATE TABLE IF NOT EXISTS `pathways` (
  `path_id` int NOT NULL AUTO_INCREMENT,
  `pathway_name` varchar(255) NOT NULL DEFAULT '',
  `info_url` varchar(255) DEFAULT NULL,
  `description_brief` varchar(255) DEFAULT NULL,
  `contributor` varchar(255) DEFAULT NULL,
  `name_systematic` varchar(255) DEFAULT NULL,
  `collection` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`path_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;

-- Insert a few sample pathways for testing
-- These should be replaced with real pathway data in production
INSERT INTO `pathways` (`path_id`, `pathway_name`, `info_url`, `description_brief`, `contributor`, `collection`) VALUES
(1, 'Sample Pathway 1', 'https://example.com/pathway1', 'A sample pathway for testing', 'PathScore Dev', 'test'),
(2, 'Sample Pathway 2', 'https://example.com/pathway2', 'Another sample pathway', 'PathScore Dev', 'test'),
(3, 'Sample Pathway 3', 'https://example.com/pathway3', 'Third sample pathway', 'PathScore Dev', 'test');
